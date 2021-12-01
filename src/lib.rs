use simplex_method::{SimplexError, Table};
use log::info;

fn check_integer(table: &Table) -> Option<usize> {
    for i in table.table.column(table.table.ncols() - 1).iter().enumerate() {
        if (i.1.round() - i.1).abs() >= 1E-13 {
            return Some(i.0);
        }
    }
    None
}

// ndarray lib does not provide a routine for swapping rows, columns etc
// incomplete software ig
fn last_rows_swap(table: &mut Table) {
    let table_row_cnt = table.table.nrows();
    for i in 0..table.table.ncols() {
        // have to swap manually to avoid double mutable borrow in mem::swap
        // imo this is better than wrapping in unsafe
        // yeah, rust is whack
        let tmp = table.table[[table_row_cnt - 2, i]];
        table.table[[table_row_cnt - 2, i]] = table.table[[table_row_cnt - 1, i]];
        table.table[[table_row_cnt - 1, i]] = tmp;
    }
}

fn get_table(mut table: Table, float_row: usize, is_le: bool) -> Table {
        let float_coeff = table.table[[float_row, table.table.ncols() - 1]];
        let mut constr = table.table.row(float_row).to_owned();
        if is_le {
            for i in constr.iter_mut() {
                *i *= -1f64;
            }
            *constr.last_mut().unwrap() += float_coeff.floor();
        } else {
            *constr.last_mut().unwrap() -= float_coeff.ceil();
        }
        table.table.push_row(constr.view()).expect("Somehow matrix dimensions are incompatible?");
        last_rows_swap(&mut table);
        table.base_var.insert(table.base_var.len() - 1,
            (table.base_var.len() + table.supp_var.len() - 1).to_string());
        table
}

fn choose_table(le_result: Result<Table, SimplexError>, ge_result: Result<Table, SimplexError>) -> Result<Table, SimplexError> {
    if le_result.is_ok() && ge_result.is_err() {
        info!("Only <= has solutions");
        return le_result;
    } else if le_result.is_err() && ge_result.is_ok() {
        info!("Only >= has solutions");
        return ge_result;
    } else if le_result.is_err() && ge_result.is_err() {
        info!("No solutions");
        return Err(simplex_method::SimplexError::NoSolutionsError);
    } else {
        let le_table = le_result.unwrap();
        let ge_table = ge_result.unwrap();
        let le_func = le_table.table[[le_table.table.nrows() - 1, le_table.table.ncols() - 1]];
        let ge_func = ge_table.table[[ge_table.table.nrows() - 1, ge_table.table.ncols() - 1]];
        // comparing by abs is fine: min(F) = -max(-F)
        if le_func.abs() > ge_func.abs() {
            info!("Choosing <=");
            return Ok(le_table);
        } else {
            info!("Choosing >=");
            return Ok(ge_table);
        }
    }
}

pub fn bnb_optimise(mut table: Table) -> Result<Table, SimplexError> {
    table.optimise()?;
    info!("Simplex optimised:\n{}", table);
    if let Some(float) = check_integer(&table) {
        info!("Found float at row {}", float);
        info!("Making <= table");
        let le_table = get_table(table.clone(), float, true);
        info!("<= table constructed:\n{}", le_table);
        let le_result = bnb_optimise(le_table);

        info!("Making >= table");
        let ge_table = get_table(table.clone(), float, false);
        info!(">= table constructed:\n{}", ge_table);
        let ge_result = bnb_optimise(ge_table);
        return choose_table(le_result, ge_result);
    }
    Ok(table)
}
