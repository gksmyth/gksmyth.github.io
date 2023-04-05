# reverse order of publication in html page

x <- readLines("2003to2014-short.html")
i <- which(x=="")
n <- length(x)
npubs <- length(i)+1
first <- c(1,i+1)
last <- c(i-1,n)

con <- file("2003to2014-short-reversed.html", "w")  
for (pub in npubs:1) {
  writeLines(x[first[pub]:last[pub]],con)
  writeLines("",con)
}

close(con)
