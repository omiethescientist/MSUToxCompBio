#R for Data Science Studying

#1/10/2019
#Loading Libraries
library(tidyverse)
mpg
#Creating a ggplot for miles per gallon data (mpg)
# displ = car engine size
# hwy = highway mpg 
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))

#So we start with  ggplot()
#Then we add layers
#geom_point for scatter plots
#geom takes a mapping arguement 
#that is paired with aes which specify the x and y axes
#so lets make a generic code
#ggplot(data = <data>) + 
#  <GEOM_FUNCTION>(mapping = aes(<Mappings>))
# 3.2.4 excerscises
ggplot(data = mpg)
length(mpg$displ)
length(colnames(mpg))
?mpg
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = cyl, y = hwy))
ggplot(data = mpg) + 
         geom_point(mapping = aes(x = drv, y = class))

#3.3 Aesthetic Mappings
#Colors and Shapes for your points to add a categorical dimension to your graph
#e.g.
ggplot(data= mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = class))
#Since the data is categorical it is important to use a categorical identifier rather
#than a quantitative one like size or a continous number line
#e.g. size
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, size = class))
#aes  = aesthetic
#ggplot can only use six shapes at a time by default

#3.3.1 excercises
#1. What is wrong with this code?
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy), color = "blue")
#Nothing, they have updated the code since writing this book

#Which variables in mpg are categorical? Which variables are continuous? (Hint: type ?mpg to read the documentation for the dataset). How can you see this information when you run mpg?
?mpg
#model, trans, drv, fl, ckass
#displ, year, cyl, city, hwy

#Map a continuous variable to color, size, and shape. How do these aesthetics behave differently for categorical vs. continuous variables?
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = cty))
#Uses color gradient
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, size = cty))
#Size Gradient
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, shape = cty))
#Error: A continuous variable can not be mapped to shape

#What happens if you map the same variable to multiple aesthetics?
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, size = cty, color = cty))
#It maps them to both aesthetics

#What does the stroke aesthetic do? What shapes does it work with? (Hint: use ?geom_point)
?geom_point
vignette("ggplot2-specs")
#Control the size of the stroke of a line (maybe?)

#What happens if you map an aesthetic to something other than a variable name, like aes(colour = displ < 5)? Note, you’ll also need to specify x and y.
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy, size = cty, color = displ < 5))
#It labels based on the boolean True/False value

#3.4 Your code won't work sometimes. It's ok.Computers are picky.

#3.5 Facets
ggplot(data = mpg) + 
  geom_point(mapping = aes(x =displ, y = hwy)) +
  facet_wrap(~ class, nrow = 2)
#facet_wrap for one variable
#facet_grid for two
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(drv ~ cyl)

ggplot(data = mpg) + 
  geom_point(mapping = aes(x =displ, y = hwy)) +
  facet_wrap(~ cty, nrow = 2)
#it tries to pretend its categorical until it can't b/c of memory

#there ar no cars with with rear wheel drive and 4 or 5 cylinders

ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy)) +
  facet_grid(.~cyl)
#it sets the axis on x or y for the single row that it assumes

#two plots
ggplot(data = mpg) + 
  geom_point(mapping = aes(x = displ, y = hwy))
ggplot(data = mpg) +
  geom_smooth(mapping = aes(x = displ, y = hwy))

#linetype aesthetic
ggplot(data = mpg) + 
  geom_smooth(mapping = aes(x = displ, y = hwy, linetype = drv))

#Overlaying plots
ggplot(data = mpg) + 
  geom_smooth(mapping = aes(x = displ, y = hwy, linetype = drv, color = drv)) + 
  geom_point(mapping = aes(x = displ, y = hwy, color = drv))

#se parameter allows you to add confidence intervals to the code

#3.7 Statistical transformations
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x = cut))

#Overiding default stat
demo <- tribble(
  ~cut, ~freq,
  "fair", 1610,
  "Good", 4906,
  "Very Good", 12082,
  "Premium", 13791,
  "Ideal", 21551
)

ggplot(data = demo) +
  geom_bar(mapping = aes(x = cut, y = freq), stat = 'identity')

#Normalized Bar Chart
ggplot(data = diamonds)+
  geom_bar(mapping = aes(x = cut, y=..prop..,group = 1))

#Stat Summary
ggplot(data = diamonds) + 
  stat_summary(
    mapping = aes(x = cut, y = depth),
    fun.ymin = min,
    fun.ymax = max,
    fun.y = median
  )
#Exercises
#1
#geom_errorbar


#3.8 Position adjustments
ggplot(data = diamonds) +
  geom_bar(mapping = aes(x=cut, fill=cut))
#aading another variable
ggplot(data = diamonds) + 
  geom_bar(mapping = aes(x=cut, fill = clarity))

#Summary of Chapter 3
#ggplot(data = <DATA>) + 
#<GEOM_FUNCTION>(
#  mapping = aes(<MAPPINGS>),
#  stat = <STAT>, 
#  position = <POSITION>
#) +
#  <COORDINATE_FUNCTION> +
#  <FACET_FUNCTION>

#Workflow: Basics
# Data transformation
library(nycflights13)
library
#When dealing with irrational numbers use near() instead of an == boolean
#filter(flights, month == 11| month == 12) is alias of filter(flights, month %in% c(11,12))