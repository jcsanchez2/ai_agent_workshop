#!/usr/bin/env Rscript

bottles <- function(n) {
  if (n == 0) "no more bottles" else if (n == 1) "1 bottle" else paste(n, "bottles")
}

for (n in 99:1) {
  cat(sprintf("%s of beer on the wall, %s of beer.\n", bottles(n), bottles(n)))
  cat(sprintf(
    "Take %s down and pass it around, %s of beer on the wall.\n\n",
    if (n == 1) "it" else "one",
    bottles(n - 1)
  ))
}

cat("No more bottles of beer on the wall, no more bottles of beer.\n")
cat("Go to the store and buy some more, 99 bottles of beer on the wall.\n")
