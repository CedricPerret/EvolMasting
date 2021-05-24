source("C:/Users/cedri/OneDrive/Research/B1-Codes/R_template/Toolbox.R")

f_plot_math("1/(1+exp(-0.05*(x-100)))",0,4000)

f_plot_math("1-1/(1+exp(-10*(x-1)))",0,10)

f_plot_math("x^6",0,1)

#For proportion of dead individuals
f_plot_math("1/(1+exp(-0.8*(x-3)))",0,10)+ylim(0,1)


wd="C:/Users/cedri/OneDrive/Research/A1-Projects/2021_EvolMasting/Code"
setwd(wd)
ls=list.files(pattern="csv")
ls
data_raw = fread(ls[2])
data_raw


ggplot(data_raw[year<50&population==1],aes(x=year,y=stock,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==1],aes(x=year,y=resources,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==1],aes(x=year,y=n_flowers,colour=ID))+
  geom_point()+
  geom_line()


ggplot(data_raw[year<50&population==5],aes(x=year,y=stock,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==5],aes(x=year,y=stock + resources,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==5],aes(x=year,y=alpha,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==5],aes(x=year,y=n_flowers,colour=ID))+
  geom_point()+
  geom_line()


ggplot(data_raw[year<50&population==2],aes(x=year,y=stock,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==2],aes(x=year,y=alpha,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50&population==2],aes(x=year,y=(stock + resources),colour=ID))+
  geom_point()+
  geom_line()


ggplot(data_raw[year<50&population==2],aes(x=year,y=n_flowers,colour=ID))+
  geom_point()+
  geom_line()

ggplot(data_raw[year<50],aes(x=population,y=n_seeds, group = population))+
  geom_boxplot()


data_raw = f_import_data(c("detail=0","beta=1","init=matching"),type_file = ".csv")
data_raw
ggplot(data_raw,aes(x=year,y=n_ind,fill=strategy))+
  geom_col(width=1)
