context("Export of results as CSV file")

test_that("CSV file structure is expected for single evidence", {
      
      res = ActivePathways(dat[,1, drop = F], gmt)
	CSV_fname = "res.csv"
      suppressWarnings(file.remove(CSV_fname))
    
	export_as_CSV(res, CSV_fname)
      res_from_CSV = read.csv(CSV_fname, stringsAsFactors = F)
      suppressWarnings(file.remove(CSV_fname))
    
      expect_equal(colnames(res_from_CSV), colnames(res))
      expect_equal(res_from_CSV$term_id, res$term_id)

})


test_that("CSV file structure is expected for multiple evidence", {

	res = ActivePathways(dat, gmt)
	CSV_fname = "res.csv"
      suppressWarnings(file.remove(CSV_fname))
    
	export_as_CSV(res, CSV_fname)
      res_from_CSV = read.csv(CSV_fname, stringsAsFactors = F)
      suppressWarnings(file.remove(CSV_fname))
    
      expect_equal(colnames(res_from_CSV), colnames(res))
      expect_equal(res_from_CSV$term_id, res$term_id)

})


test_that("CSV file is exported when there are NULL entries", {
      
      res = ActivePathways(scores=scores_test, gmt=gmt_reac, significant=1)
      CSV_fname = "res.csv"
      suppressWarnings(file.remove(CSV_fname))
      
      export_as_CSV(res, CSV_fname)
      res_from_CSV = read.csv(CSV_fname, stringsAsFactors = F)
      suppressWarnings(file.remove(CSV_fname))
      
      # convert res overlap column values to a string type
      res$overlap <- sapply(res$overlap, function(x) paste(x, collapse = "|"))
      expect_equal(res_from_CSV$overlap, res$overlap)

})
