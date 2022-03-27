#install.packages("tidyverse")

library(tidyverse)
detach("package:tidyverse", unload = TRUE)
str_c("hello")
#require(X)
help(tidyverse) # but also googling, stackoverflow.
?tidyverse
?rep
# Part 1: Object types & class ####
### create & view an object----
my_value <- 24
(my_value <- 24)
my_value <- 24;my_value
my_value
my_vector
my_string <- "24"

my_mat <- matrix(data=1:9,nrow=3)
length(my_mat)
dim(my_mat)
str(my_mat)
summary(my_mat)

typeof(my_mat)
class(my_mat)
typeof(my_value)
class(my_value)
typeof(my_string)
class(my_string)

rm(list = ls())

# c function (c stands for "combine")
my_vector <- c(5,34,76,13)
my_vector <- 35:50
my_vector = 1:10
1:10
rep(1,5,2) #?rep
rep(1,100,2)
# Exercise: create a vector
my_vector[5]


## Removing Variables, Setting Working Directory ##


# Set your working directory -- This is where R goes to look for files and save
# stuff by default. You will need to do this for each computer you run your
# script file on. In RStudio, you can go to 'Session' -> 'Set Working Directory'
# -> 'Choose Directory', and select a folder from a drop down menu. For me, this
# looks like:
setwd("~/Desktop")



# Part 2: Comparison & subsetting ------------------------

## The Basic Operators ##

5 < 6 # 5 is less than 6: returns TRUE
5 > 6 # 5 is not greater than 6: returns FALSE
5 == 5 # 5 is equal to 5: returns TRUE
5 != 6 # 5 is not equal to 6: returns TRUE
5 <= 5 # 5 is less than or equal to 5: returns TRUE

# We can index the elements of matrices in a variety of ways:
iris
iris[1,] # rows come before the comma [row,colum]
iris[,5] # columns come after
iris[3,3]

# We can create new variables out of pieces of matrices:
val <- my_matrix[3,3]
val <- my_matrix[,3]


## Exercise: How do we achieve the same goal with the tidyverse package? ####

# Part 3: Data import ####
write.csv(iris,"C:/Users/xwuey/Desktop/iris.csv")

df <- read.csv("C:/Users/xwuey/Desktop/iris.csv") 
df <- tibble(df)
head(df,20)
class(df)
data.frame(names(df))
names(df) %>% data.frame()
nrow(df)
ncol(df)
dim(df)
str(df)
summary(df)

# Subsetting columns (var's) & rows----
## base R
df[, "Sepal.Length"]
df[,c("Sepal.Length","Sepal.Width")]
df$Sepal.Length
df[["Sepal.Length"]]
df[2]
df[c(1, 2, 3), ] #get first 3 rows
df[c(1:20),]
df[which(df$Sepal.Length>1 & df$Sepal.Width>2),] # get rows by criteria
df$Petal.Width

## tidyverse
df %>% filter(Sepal.Length>1,Sepal.Width>2)
filter(df,Sepal.Length==5.1|Sepal.Length==4.4)
filter(df,Sepal.Length %in% c(5.1,4.4))

a <- df %>% filter(Sepal.Length>1,Sepal.Width>2) %>% 
  filter(X>5) #%>% 
  
# get col by criteria
df %>% select(Sepal.Length,Sepal.Width)

# identical results?
df[which(df$Sepal.Length>1),]==df %>% filter(Sepal.Length>1) 
identical(df[which(df$Sepal.Length>1),],df %>% filter(Sepal.Length>1))
sum(df[which(df$Sepal.Length>1),]!=df %>% filter(Sepal.Length>1)) # count using the property of T=1/F=0

# Exercise: show the values of Species & Petal.Width for flowers with Sepal.Length larger or equal to 3, using both base R and tidyverse grammars

# mutate and summarise
df$new_var = df$Sepal.Length*10
df2 <- df %>% dplyr::mutate(new_var2 = c(1:150))
mean(df$Sepal.Length)
summary(df[,c(2:4)])
a <- df %>% group_by(Species) %>% mutate(m = mean(Sepal.Length))
a <- df %>% 
  group_by(Species) %>% 
  summarise(m = mean(Sepal.Length),
            k = median(Sepal.Length)) %>% 
  ungroup()
