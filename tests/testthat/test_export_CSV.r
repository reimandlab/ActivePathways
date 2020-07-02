context("Export of results as CSV file")

test_that("CSV file structure is expected for single evidence", {

	res = ActivePathways(dat[,1, drop = F], gmt)
	CSV_fname = "res.csv"
    suppressWarnings(file.remove(CSV_fname))
    
	export_as_CSV(res, CSV_fname)
    res_from_CSV = read.csv(CSV_fname, stringsAsFactors = F)
    suppressWarnings(file.remove(CSV_fname))
    
    expect_equal(colnames(res_from_CSV), colnames(res))
    expect_equal(res_from_CSV$term.id, res$term.id)

})


test_that("CSV file structure is expected for multiple evidence", {

	res = ActivePathways(dat, gmt)
	CSV_fname = "res.csv"
    suppressWarnings(file.remove(CSV_fname))
    
	export_as_CSV(res, CSV_fname)
    res_from_CSV = read.csv(CSV_fname, stringsAsFactors = F)
    suppressWarnings(file.remove(CSV_fname))
    
    expect_equal(colnames(res_from_CSV), colnames(res))
    expect_equal(res_from_CSV$term.id, res$term.id)

})