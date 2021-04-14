source("C:/Users/cedri/OneDrive/Research/B1-Codes/R_template/Toolbox.R")

f_plot_math("1/(1+exp(-0.05*(x-100)))",0,4000)

f_plot_math("1/(1+exp(-0.003*(x-2000)))",0,4000)

wd="C:/Users/cedri/OneDrive/Research/B1-Codes"
setwd(wd)
list.files(pattern="csv")
data_raw = fread("valentin_trees_res.csv")

data_raw
pt_Y=ggplot(data_raw[individual<5],aes(x=gen,y=Y, group = individual,colour=as.factor(individual)))+
  geom_point()+
  geom_line()+
  u_theme1
pt_Y

pt_pollinisation=ggplot(data_raw[individual<5],aes(x=gen,y=pollination_rate, group = individual,colour=as.factor(individual)))+
  geom_point()+
  geom_line()+
  u_theme1
pt_pollinisation

plot_grid(pt_Y,pt_pollinisation,nrow=2)
