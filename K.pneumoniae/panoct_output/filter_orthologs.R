library(dplyr)
library(tidyverse)
library(ggplot2)

# Read in match table
match.table <- read.delim(header = FALSE, "matchtable.txt")
colnames(match.table) <- c("ortholog", paste("strain", 1:74, sep=""))

# Determine number of genes per ortholog
counter = 1
num.strains <- c()
while(counter <= dim(match.table)[1])
{
  col.num = 2
  valid.num = 0
  while(col.num <= 75)
  {
    if(as.vector(match.table[counter,col.num]) != "----------")
    {
      valid.num = valid.num+1
    }
    col.num <- col.num+1
  }
 num.strains <- append(num.strains, valid.num)
 counter <- counter+1
}
match.table <- cbind(match.table, num.strains)
names(match.table) [76] <- "total"

# Filter poor matches
hist(match.table$total)
cutoff = 0
match.table %>% 
  filter(total > cutoff) %>% 
  summarise(sum(total)) %>% 
  ggplot(.,aes(x=total))+geom_line(stat="density")


cutoff.vect <- c()
sum.vect <- c()
while(cutoff < 75)
{
  sum.vect <- append(sum.vect,match.table %>% 
    filter(total > cutoff) %>% 
    summarise(sum(total)))
  cutoff.vect <- append(cutoff.vect, cutoff)
  cutoff <- cutoff +10
}
plot(cutoff.vect, sum.vect)









