db_w_microarray <- 
function(
	df, 
	dfname,
	dbname = "microarraydb"	
){
	library(RSQLite) # read library
	driver = dbDriver("SQLite") # driver setup
	db_wd <- paste("C:/SQLite/", dbname, sep="")
	con = dbConnect(driver, db_wd) # Connection object setup

	rs=dbWriteTable(con, dfname, df, row.names=T) # write dataframe as table to database
}

db_r_microarray <-
function(
	dfname,
	dbname = "microarraydb"
){
	library(RSQLite) # read library
	db_wd <- paste("C:/SQLite/", dbname, sep="")
	con = dbConnect(driver, db_wd) # Connection object setup
	temp <- paste("SELECT * FROM", dfname, sep=" ")
	return(dbGetQuery(con, temp))
	}
	
