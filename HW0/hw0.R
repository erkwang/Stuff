#STA 250 HOMEWORK 0
#Yichuan Wang

#Problem 1
count = rep(0, 100)
for (i in 1:100){
  if (!(as.integer(i/3)-i/3)) count[i] = count[i]+1
  if (!(as.integer(i/5)-i/5)) count[i] = count[i]+2
}
foo = as.character(1:100)
foo[count == 1] = "Fizz"
foo[count == 2] = "Buzz"
foo[count == 3] = "FizzBuzz"
cat(foo)

#Problem 2
x = runif(n = 100, min=0, max=2*pi)
y = runif(n = 100, min=0, max=1)
pairs = data.frame(x = rep(x, 100), y = rep(y, each = 100))
pairs$u = with(pairs, y*cos(x))
pairs$v = with(pairs, y*sin(x))
plot(v~u, data = pairs)
pairs$r = with(pairs, sqrt(u^2+v^2))

#Problem 3
snippet = "Hello, my name is Bob. I am a statistician. I like statistics very much."
snippet = strsplit(snippet, split = "")[[1]]
for (i in 1:length(snippet)){
  writeLines(text = snippet[i], sep = "",
             con = paste("~/Desktop/SkyDrive/STA 250/hw0/out_", i, ".txt", sep = ""))
}
outfiles = list.files("~/Desktop/SkyDrive/STA 250/hw0/", pattern = "^out_", full.names=TRUE)
library(gtools)
outfiles = mixedsort(outfiles)
snipout = sapply(outfiles, function(f)readLines(con = f, n = -1))
snipout = paste(snipout, collapse = "")
writeLines(snipout, sep = "", con = ("SkyDrive/STA 250/hw0/snippet.txt"))

#Prolem 6









