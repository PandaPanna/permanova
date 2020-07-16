d <- read.csv('d.csv')
dim(d)
# There are 307 samples, 12387 otus
library(tidyr)


# Check if there is missing columns

# Method 1: check if 'Otu8' is in col_names
col_names <- colnames(d)
'Otu8' %in% col_names


# Method 2: sort columns in numeric order
col_index <- c()
for (val in col_names)
{
  # Slice the column names from the 4th letter
  # e.g. 'Otu8' will be '8'
  col_index <- c(col_index, substr(val,4, 100))
}
sort(as.numeric(col_index))


# Change data to a dataframe
df<- data.frame(d)
head(df)
# Transform the dataframe to create the frequency table
total_df <- df %>% gather(Otu, Frequency, col_names[2:12388])
head(total_df,10)
colnames(total_df)[1:2] <- c('Var1', 'Var2')
dim(total_df)# 3802809 rows
12387*307# Match


# Check the otu8 
check_otu8 <- total_df[total_df[, "Var2"] == 'Otu8',]
head(check_otu8)

# How many unique Otus have occured in total?
length(unique(total_df[,'Var2']))
# Matches the total Otus in the database 12387

# How many unique Otus have occured in Blank?
bl_names <- c('BL1', 'BL2','BL3','BL4','BL5')
# If you want to include negative samples just create a name list including BL and Neg...

# We only need the data match the condition: Var1 value in the bl_names list.
bl_df <- total_df[total_df$Var1 %in% bl_names & total_df$Frequency !=0,]
length(unique(bl_df[,'Var2']))
# There are 3421 unique Otus being found in blank samples.
3421/12387
# 27.62 % of the total Otus can be found in blank.
dim(bl_df)
# 3421 Otus have occured 5502 times in 5 blank samples.

head(total_df)
# We need another column which is the 'Max'

# I transposed the dataframe first in permanova_04, it's easier
# to find the max value in a row first, then spread the dataframe into frequency table.
# Here since we already got the frequency table total_df, 
# we will have to gather the same 'Var2' value to find the max for each Otu group


# First, calculate the max value by each Otu group.(One unique Var2 has one max value) 
Max <- tapply(total_df$Frequency, total_df$Var2,max)
head(Max)
length(Max)

# Second, assign max value to each row.
total_df <- data.frame(total_df, max=rep(Max, table(total_df$Var2)))
head(total_df,100)
total_df[total_df$Var2=='Otu1',]


######Linear Model######
# Assumption1: as max increases, the frequency in bl should also increase.
# Now the total_df table is still 3802809 rows
# If we use the total table to build linear model:

lm_total <- lm(total_df$Frequency~total_df$max)
summary(lm_total)
# You can see that the intercept and slope are all 0 which is meaningless
# It's because most of the Otus never occurs in Blank Samples.

# Thus we have to get rid of those otus which never occur in BL
bl_names <- c('BL1', 'BL2','BL3','BL4','BL5')
y_threshold = 5 # If an Otu occurence in BL is less than 6, we will consider it's zero. 
bl_df <- total_df[total_df$Var1 %in% bl_names & total_df$Frequency > y_threshold,]



head(bl_df) 
dim(bl_df)
# 2597 points now
# As we mentioned in line 49
# 3421 Otus have occured 5502 times in blank samples.
# Since some of the occurence is less than 6, so the dataset shrinks after 
# using shreshold =5 as the upper limit

# The top 5  occurences in BL
bl_df[order(bl_df$Frequency, decreasing = TRUE),][c(1:5),]
# The least 5  occurences in BL
bl_df[order(bl_df$Frequency, decreasing = FALSE),][c(1:5),]

#Top 10 max value
bl_df[order(bl_df$max, decreasing = TRUE),][c(1:10),]
#Otu332 has has max abundance 59174 which is very high,
#but it only has abundance in BL 28~89 (89/59174=0.0015), I think we can ingnore it
#right??????
plot(bl_df$max, bl_df$Frequency, col=ifelse(bl_df$Var2=='Otu332', "red", "black"))


# Method1: you can limit the x asix threshold to get rid of thoese red points.
# This might not be the best idea...
x_threshold <- 10000
plot(bl_df$max, bl_df$Frequency, col=ifelse(bl_df$max>x_threshold, "red", "black"))

# Method2: if y/x <0.01, get rid off the points.
plot(bl_df$max, bl_df$Frequency, col=ifelse(bl_df$Frequency/bl_df$max<0.01, "red", "black"))


bl_df <- bl_df[bl_df$Frequency/bl_df$max >= 0.01,]
dim(bl_df) #2497 observations
length(unique(bl_df$Var2)) #1655 Otus to be considered


plot(bl_df$max, bl_df$Frequency)
bl_lm <- lm(bl_df$Frequency~bl_df$max)
summary(bl_lm)
abline(bl_lm, col= 'red')



###### Explore: top 5 most occured otus in bl###### 
# Add bl_max column
bl_max <- tapply(bl_df$Frequency, bl_df$Var2,max)
bl_df <- data.frame(bl_df, bl_max=rep(bl_max, table(bl_df$Var2)))
top <- bl_df[order(bl_df$bl_max, decreasing = TRUE),]
top <- top[!duplicated(top$Var2),][c(1:5),]
top5_otus <- c(top['Var2'])
top5_otus

ggplot(data=bl_df,aes(x=max,y=Frequency))+
  geom_point(color='#00AFBB') +
  geom_point(data=bl_df[bl_df$Var2=='Otu8',],color="#FC4E07",size=2)+
  geom_point(data=bl_df[bl_df$Var2=='Otu6',],color="red",size=2)+
  geom_point(data=bl_df[bl_df$Var2=='Otu22',],color="yellow",size=2)+
  geom_point(data=bl_df[bl_df$Var2=='Otu9',],color="green",size=2)+
  geom_point(data=bl_df[bl_df$Var2=='Otu35',],color="purple",size=2)+
  geom_abline(slope=bl_lm$coefficients[2],intercept=bl_lm$coefficients[1], colour = '#E7B800')+
  geom_quantile(quantiles = 0.95,colour = "purple", size = 2, alpha = 0.5)





######Plot beautifully  ######
library(ggplot2)
library(quantreg)

# The color code in R is:
# http://www.sthda.com/english/wiki/colors-in-r
# You can play with it to change the color as you like
# Both color="#FF0000" and color="red" are totally the same thing
# They are the same red color
# You don't have to dig too deep into the color 
# I use #FC4E07 #00AFBB etc just because they look more beautiful in the graph than default color



graph1 <- ggplot(data=bl_df,aes(x=max,y=Frequency))+
  #Plot all the data point with color'#00AFBB' which is a baby blue color.
  geom_point(color='#00AFBB', size=1) +
  # Show the Otu8 data with different color '#FC4E07'
  geom_point(data=bl_df[bl_df$Var2=='Otu8',],color="#FC4E07",size=1)+
  # Add text label(Var1: sample names) to Otu8 observations
  geom_text(data=subset(bl_df, Var2=='Otu8'),
            aes(label=Var1))+
  # Add bl_lm line as yellow color to the plot
  geom_abline(slope=bl_lm$coefficients[2],intercept=bl_lm$coefficients[1], colour = '#E7B800')

graph1


#Quantile regression
bl_rq <- rq(bl_df$Frequency~bl_df$max,tau = 0.95)
bl_rq
# Get the coefficients so we can filter out all unwanted points 

graph1 <- graph1+
  # Add bl_lm line as red color to the plot
  geom_abline(slope=bl_rq$coefficients[2],intercept=bl_rq$coefficients[1], colour = 'red')

graph1
# To see the data above bl_rq line
# Those point match the codition y > 0.366667*x + 16.46667 will be painted red
graph1 <- graph1+geom_point(data=bl_df[bl_df$Frequency> 0.366667*bl_df$max + 16.46667,],color="red",size=1)
graph1


# Show some Otu name label when y > 450
graph1 <- graph1+geom_text(data=subset(bl_df, (bl_df$Frequency> 0.366667*bl_df$max + 16.46667)& bl_df$Frequency>450),
                    aes(label=Var2))
graph1

# above_obs is all the points above the bl_rq line  y > 0.366667*x + 16.46667
above_obs <- bl_df$Var2[bl_df$Frequency> (0.366667*bl_df$max+ 16.46667)]
length(above_obs)
# 120 observations are above the bl_rq line

abnormal_otus <- unique(above_obs)
length(abnormal_otus)
abnormal_otus
# There are 111 otus to be considered as abnormal,


# Color all the abnormal otus, both above and under the bl_rq line
# In the bl_df dataframe, If the Var2 value is in the abnormal_otus list, then it wil be colored "#FC4E07"
graph1 <- graph1+geom_point(data=bl_df[bl_df$Var2 %in% abnormal_otus,],color="#FC4E07",size=1)
graph1

# Zoom in to see detail
graph1 + xlim(0, 2000) + ylim(0,500)

