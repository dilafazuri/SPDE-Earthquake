#==================== STEPS SPDE ========================================================================================
#========================================================================================================================
#==================== 1. Install and Read the Library
#========================================================================================================================
setwd("/Users/dila.azuri/Documents/Magister/Tesis/SPDE/Syntax")
library("devtools")
library(ggplot2)
library(sf)
library(gstat)
library(sp) 
library(INLA)
library(writexl)
library(dplyr)
library(plotly)
library(akima)
library(cowplot)
library(patchwork)

#========================================================================================================================
#==================== 2. Input Data
#========================================================================================================================
#===========2.1 Input Data Earthquake
Earthquakes<-read.csv("INLA - Gempa.csv",sep=";")
head(Earthquakes)
Earthquakes20=Earthquakes[Earthquakes$Tahun==2020,]
head(Earthquakes20)
Earthquakes20$LAT <- gsub(",", ".", Earthquakes20$LAT)
Earthquakes20$LAT <- as.numeric(Earthquakes20$LAT)
Earthquakes20$LON <- gsub(",", ".", Earthquakes20$LON)
Earthquakes20$LON <- as.numeric(Earthquakes20$LON)
Earthquakes20$MAG <- gsub(",", ".", Earthquakes20$MAG)
Earthquakes20$MAG <- as.numeric(Earthquakes20$MAG)
colnames(Earthquakes20)[colnames(Earthquakes20) == "Bulan"] <- "Month"
head(Earthquakes20)

#=========2.1.1 Exploratory Data Analysis
#========2.1.1.1 Statistics Descriptive
mag_stats_by_month <- Earthquakes20 %>%
  group_by(Month) %>%
  summarise(
    count = n(),
    mean_mag = mean(MAG, na.rm = TRUE),
    median_mag = median(MAG, na.rm = TRUE),
    sd_mag = sd(MAG, na.rm = TRUE),
    min_mag = min(MAG, na.rm = TRUE),
    max_mag = max(MAG, na.rm = TRUE)
  )
mag_stats_by_month

mean(Earthquakes20$MAG)
sd(Earthquakes20$MAG)
min(Earthquakes20$MAG)
max(Earthquakes20$MAG)

#========2.1.1.2 Boxplot
ggplot(Earthquakes20, aes(x = as.factor(Month), y = MAG, fill = as.factor(Month))) +   
  geom_boxplot(color = "black") +   
  labs(title = "Boxplot of Earthquake Magnitude by Month",
       x = "Month",
       y = "Magnitude (Mw)") +   
  scale_x_discrete(labels = c("Jan", "Feb", "Mar", "Apr", "May", "June", 
                              "July", "Aug", "Sept", "Oct", "Nov", "Dec")) +
  scale_fill_manual(values = c("1" = "cadetblue3", "2" = "#377EB8", "3" = "cornsilk1",
                               "4" = "#984EA3", "5" = "coral", "6" = "deeppink3",
                               "7" = "bisque", "8" = "#F781BF", "9" = "cornsilk4",
                               "10" = "aliceblue", "11" = "brown", "12" = "darkseagreen")) +
  theme_minimal()+
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    panel.border = element_rect(color = "white", fill = NA, linewidth = 1),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid")
  )

#===========2.2 Input Sumatra Maps
Indonesia<-readRDS('gadm36_IDN_1_sp.rds') 
Indonesia$NAME_1
Sumatra<-Indonesia[Indonesia$NAME_1 %in%c("Kepulauan Riau","Sumatera Selatan","Lampung","Sumatera Utara","Bangka Belitung","Jambi","Riau","Bengkulu","Sumatera Barat","Aceh"),]
Sumatra_sf <- st_as_sf(Sumatra)

#===========2.3 Input Active Faults
fault = st_read(file.choose())
sumatera_fault <- st_bbox(c(xmin = 95, ymin = -6, xmax = 106, ymax = 6), crs = st_crs(fault))
fault_sumatera <- st_crop(fault, sumatera_fault)

#===========2.4 Epicenters Visualization
Month_names <- c( "1"="January",
                  "2"="February",
                  "3"="March",
                  "4"="April",
                  "5"="May",
                  "6"="June",
                  "7"="July",
                  "8"="August",
                  "9"="September",
                  "10"="October",
                  "11"="November",
                  "12"="December")   

ggplot() +
  geom_sf(data = Sumatra_sf) +
  geom_sf(data = fault_sumatera, color = "blue", size = 1) +
  geom_point(data = as.data.frame(Earthquakes20), aes(x = LON, y = LAT, color =MAG)) +  
  scale_color_gradient(low = "yellow", high = "red")+
  labs(title = "Distribution of relocated earthquake epicenters", x = "Longitude", y = "Latitude") +
  theme_bw()

#========================================================================================================================
#==================== 3. Data Preparation
#========================================================================================================================
locations <- unique(Earthquakes20[, c("LON", "LAT")])
all_combinations <- merge(locations, data.frame(Month = unique(Earthquakes20$Month)))
head(all_combinations)
final_Earthquakes20 <- merge(all_combinations, Earthquakes20, 
                             by = c("LON", "LAT", "Month"), all.x = TRUE)
head(final_Earthquakes20)

final_Earthquakes20$MAG[is.na(final_Earthquakes20$MAG)] <- NA
final_Earthquakes20$DEPTH[is.na(final_Earthquakes20$DEPTH)] <- NA
final_Earthquakes20$Year[is.na(final_Earthquakes20$Year)] <- 2020

Earthquakes20<-data.frame(LON=final_Earthquakes20$LON,LAT=final_Earthquakes20$LAT, 
                          Month=final_Earthquakes20$Month, MAG=final_Earthquakes20$MAG, 
                          DEPTH=final_Earthquakes20$DEPTH)
head(Earthquakes20)

#========================================================================================================================
#==================== 4. Grid
#========================================================================================================================
#===========4.1 Several Interval Grid
polygon_coords <- data.frame(
  longitude = c(94, 94, 98, 108, 102, 94),
  latitude = c(2, 6, 6, -7, -7, 2)
)

polygon_sf <- st_as_sf(polygon_coords, coords = c("longitude", "latitude"), crs = 4326)
polygon_utm <- st_transform(polygon_sf, crs = 32648)
utm_coords <- as.data.frame(st_coordinates(polygon_utm))
utm_coords_df <- data.frame(longitude = utm_coords$X, latitude = utm_coords$Y)

#=========4.1.1 10km x 10km
interval1 <- 10000
grid_points1 <- expand.grid(
  longitude = seq(min(utm_coords$X), max(utm_coords$X), by = interval1),
  latitude = seq(min(utm_coords$Y), max(utm_coords$Y), by = interval1)
)

inside1 <- point.in.polygon(grid_points1$longitude, grid_points1$latitude, 
                            utm_coords$X, utm_coords$Y)
inside_points1 <- grid_points1[inside1 == 1, ]

ggplot() +
  geom_path(data = utm_coords_df, aes(x = longitude, y = latitude), color = "red") +
  geom_point(data = inside_points1, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  labs(x = "Longitude", y = "Latitude", title = "Interval Grid: 10000") +
  theme_minimal()

#=========4.1.2 20km x 20km
interval2 <- 20000
grid_points2 <- expand.grid(
  longitude = seq(min(utm_coords$X), max(utm_coords$X), by = interval2),
  latitude = seq(min(utm_coords$Y), max(utm_coords$Y), by = interval2)
)

inside2 <- point.in.polygon(grid_points2$longitude, grid_points2$latitude, 
                            utm_coords$X, utm_coords$Y)
inside_points2 <- grid_points2[inside2 == 1, ]

ggplot() +
  geom_path(data = utm_coords_df, aes(x = longitude, y = latitude), color = "red") +
  geom_point(data = inside_points2, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  labs(x = "Longitude", y = "Latitude", title = "Interval Grid: 20000") +
  theme_minimal()

#=========4.1.3 30km x 30km
interval3 <- 30000
grid_points3 <- expand.grid(
  longitude = seq(min(utm_coords$X), max(utm_coords$X), by = interval3),
  latitude = seq(min(utm_coords$Y), max(utm_coords$Y), by = interval3)
)

inside3 <- point.in.polygon(grid_points3$longitude, grid_points3$latitude, 
                            utm_coords$X, utm_coords$Y)
inside_points3 <- grid_points3[inside3 == 1, ]

ggplot() +
  geom_path(data = utm_coords_df, aes(x = longitude, y = latitude), color = "red") +
  geom_point(data = inside_points3, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  labs(x = "Longitude", y = "Latitude", title = "Interval Grid: 30000") +
  theme_minimal()

#=========4.1.4 40km x 40km
interval4 <- 40000
grid_points4 <- expand.grid(
  longitude = seq(min(utm_coords$X), max(utm_coords$X), by = interval4),
  latitude = seq(min(utm_coords$Y), max(utm_coords$Y), by = interval4)
)

inside4 <- point.in.polygon(grid_points4$longitude, grid_points4$latitude, 
                            utm_coords$X, utm_coords$Y)
inside_points4 <- grid_points4[inside4 == 1, ]

ggplot() +
  geom_path(data = utm_coords_df, aes(x = longitude, y = latitude), color = "red") +
  geom_point(data = inside_points4, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  labs(x = "Longitude", y = "Latitude", title = "Interval Grid: 40000") +
  theme_minimal()

#=========4.1.5 50km x 50km
interval5 <- 50000
grid_points5 <- expand.grid(
  longitude = seq(min(utm_coords$X), max(utm_coords$X), by = interval5),
  latitude = seq(min(utm_coords$Y), max(utm_coords$Y), by = interval5)
)

inside5 <- point.in.polygon(grid_points5$longitude, grid_points5$latitude, 
                            utm_coords$X, utm_coords$Y)
inside_points5 <- grid_points5[inside5 == 1, ]

ggplot() +
  geom_path(data = utm_coords_df, aes(x = longitude, y = latitude), color = "red") +
  geom_point(data = inside_points5, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  labs(x = "Longitude", y = "Latitude", title = "Interval Grid: 50000") +
  theme_minimal()

#===========4.2 Grid align with the observation
#=========4.2.1 10km x 10km
coordinates(inside_points1)<-~longitude+latitude
gridded(inside_points1) = TRUE
coordsGrid.UTM1<-inside_points1
UTM<-CRS("+proj=utm +zone=48 ellps=WGS84")
Sumatra_UTM <- spTransform(Sumatra, UTM)
coordsGrid.UTM_df1 <- as.data.frame(coordsGrid.UTM1)
Sumatra_UTM_sf <- st_as_sf(Sumatra_UTM)

ggplot() +
  geom_sf(data = Sumatra_UTM_sf, fill = "white", color = "black") +
  geom_point(data = coordsGrid.UTM_df1, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  theme_minimal() +
  ggtitle("Sumatra Maps with Grid UTM")

#=========4.2.2 20km x 20km
coordinates(inside_points2)<-~longitude+latitude
gridded(inside_points2) = TRUE
coordsGrid.UTM2<-inside_points2
coordsGrid.UTM_df2 <- as.data.frame(coordsGrid.UTM2)

ggplot() +
  geom_sf(data = Sumatra_UTM_sf, fill = "white", color = "black") +
  geom_point(data = coordsGrid.UTM_df2, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  theme_minimal() +
  ggtitle("Sumatra Maps with Grid UTM")

#=========4.2.3 30km x 30km
coordinates(inside_points3)<-~longitude+latitude
gridded(inside_points3) = TRUE
coordsGrid.UTM3<-inside_points3 
coordsGrid.UTM_df3 <- as.data.frame(coordsGrid.UTM3)

ggplot() +
  geom_sf(data = Sumatra_UTM_sf, fill = "white", color = "black") +
  geom_point(data = coordsGrid.UTM_df3, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  theme_minimal() +
  ggtitle("Sumatra Maps with Grid UTM")

#=========4.2.4 40km x 40km
coordinates(inside_points4)<-~longitude+latitude
gridded(inside_points4) = TRUE
coordsGrid.UTM4<-inside_points4
coordsGrid.UTM_df4 <- as.data.frame(coordsGrid.UTM4)

ggplot() +
  geom_sf(data = Sumatra_UTM_sf, fill = "white", color = "black") +
  geom_point(data = coordsGrid.UTM_df4, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  theme_minimal() +
  ggtitle("Sumatra Maps with Grid UTM")

#=========4.2.5 50km x 50km
coordinates(inside_points5)<-~longitude+latitude
gridded(inside_points5) = TRUE
coordsGrid.UTM5<-inside_points5
coordsGrid.UTM_df5 <- as.data.frame(coordsGrid.UTM5)

ggplot() +
  geom_sf(data = Sumatra_UTM_sf, fill = "white", color = "black") +
  geom_point(data = coordsGrid.UTM_df5, aes(x = longitude, y = latitude), color = "blue", size = 0.5) +
  theme_minimal() +
  ggtitle("Sumatra Maps with Grid UTM")

#========================================================================================================================
#==================== 5. Prepare the data in the UTM system
#========================================================================================================================
st_crs(Sumatra_sf)
Sumatra_sf_utm1 <- st_transform(Sumatra_sf, crs = 32648)

Observed<-data.frame(X=Earthquakes20$LON,Y=Earthquakes20$LAT) 
coordinates(Observed) <- ~X+Y 
proj4string(Observed) <- CRS("+proj=longlat +datum=WGS84")
Observed_utm <- spTransform(Observed, CRS("+proj=utm +zone=48 +datum=WGS84 +units=m +no_defs"))
Observed_utm1<-coordinates(Observed_utm)
colnames(Observed_utm1)<-c("X","Y")

Earthquakes20$UTM_x<-Observed_utm1[,1]
Earthquakes20$UTM_y<-Observed_utm1[,2]
Earthquakes20$IT<-Earthquakes20$Month
head(Earthquakes20)

#========================================================================================================================
#==================== 6. Construct Model
#========================================================================================================================
#===========6.1 Define the Triangulation Mesh
Earthquakes_mesh <- inla.mesh.2d(
  loc=cbind(Earthquakes20$UTM_x,Earthquakes20$UTM_y), 
  offset=c(50000, 200000),
  max.edge=c(100000, 1000000))

plot(Earthquakes_mesh,col="gray80", main="Earthquakes Mesh of Sumatra")
points(Observed_utm1,col="red")
plot(Sumatra_UTM,add=T, border="blue")

#===========6.2 Define the model for observed and grids areas
coordinates.allyear<-as.matrix(data.frame(x=Earthquakes20$UTM_x,y=Earthquakes20$UTM_y))
A_est <- inla.spde.make.A(mesh=Earthquakes_mesh,loc=coordinates.allyear,
                          group=Earthquakes20$Month, n.group=12)
dim(A_est)

#===========6.3 Define the INLA parameters required for computation
control <- list(
  results = list(return.marginals.random = TRUE, return.marginals.predictor=TRUE),
  compute = list(hyperpar=TRUE, return.marginals.predictor=TRUE, return.marginals=TRUE, dic=TRUE, mlik = TRUE, cpo = TRUE, 
                 po = TRUE, waic=TRUE))

#===========6.4 Define prior distribution of Penalized Complexity
#=========6.4.1 For 10km x 10km
Grid.CoordUTM1<-as.data.frame(coordsGrid.UTM1)
m1<-nrow(Grid.CoordUTM1)
Grid.CoordUTM1<-as.matrix(Grid.CoordUTM1)
GroupTime1<-rep(c(1:12),each=m1)
All.Grid1<-kronecker(matrix(1, 12, 1), Grid.CoordUTM1)
#=========6.4.1.1 Range: 10km and Standard Deviation 0.1
prior.median.sd1001 = 0.1; prior.median.range10 = 10000

Earthquakes_spde1001 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range10, .5), 
                                           prior.sigma = c(prior.median.sd1001, .5))

s_index1001 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde1001$n.spde,
                                    n.group=12)

stack_est1001 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index1001) , tag="est")

A_pred1 <- inla.spde.make.A(mesh=Earthquakes_mesh,
                            loc=as.matrix(All.Grid1),
                            group=GroupTime1,  #selected day for prediction
                            n.group=12)

Cov1<-data.frame(UTM_x=All.Grid1[,1]/1000,UTM_y=All.Grid1[,2]/1000,DEPTH=NA,Time=GroupTime1)

stack_pred1001 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred1),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid1)),
                                                     Time=GroupTime1, UTM_x=Cov1$UTM_x, 
                                                     UTM_y=Cov1$UTM_y, DEPTH=Cov1$DEPTH), s_index1001), 
                             tag="pred")
stack1001 <- inla.stack(stack_est1001, stack_pred1001) 

#=========6.4.1.2 Range: 10km and Standard Deviation 0.5
prior.median.sd1005 = 0.5; prior.median.range10 = 10000

Earthquakes_spde1005 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range10, .5), 
                                           prior.sigma = c(prior.median.sd1005, .5))

s_index1005 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde1005$n.spde,
                                    n.group=12)

stack_est1005 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index1005) , tag="est")

stack_pred1005 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred1),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid1)),
                                                     Time=GroupTime1, UTM_x=Cov1$UTM_x, 
                                                     UTM_y=Cov1$UTM_y, DEPTH=Cov1$DEPTH), s_index1005), 
                             tag="pred")
stack1005 <- inla.stack(stack_est1005, stack_pred1005)

#=========6.4.1.3 Range: 10km and Standard Deviation 1
prior.median.sd101 = 1; prior.median.range10 = 10000

Earthquakes_spde101 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range10, .5), 
                                          prior.sigma = c(prior.median.sd101, .5))

s_index101 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde101$n.spde,
                                   n.group=12)

stack_est101 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index101) , tag="est")

stack_pred101 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred1),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid1)),
                                                    Time=GroupTime1, UTM_x=Cov1$UTM_x, 
                                                    UTM_y=Cov1$UTM_y, DEPTH=Cov1$DEPTH), s_index101), 
                            tag="pred")
stack101 <- inla.stack(stack_est101, stack_pred101) 

#=========6.4.1.4 Range: 10km and Standard Deviation 1.5
prior.median.sd1015 = 1.5; prior.median.range10 = 10000

Earthquakes_spde1015 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range10, .5), 
                                           prior.sigma = c(prior.median.sd1015, .5))

s_index1015 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde1015$n.spde,
                                    n.group=12)

stack_est1015 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index1015) , tag="est")

stack_pred1015 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred1),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid1)),
                                                     Time=GroupTime1, UTM_x=Cov1$UTM_x, 
                                                     UTM_y=Cov1$UTM_y, DEPTH=Cov1$DEPTH), s_index1015), 
                             tag="pred")
stack1015 <- inla.stack(stack_est1015, stack_pred1015) 

#=========6.4.1.5 Range: 10km and Standard Deviation 2
prior.median.sd102 = 2; prior.median.range10 = 10000

Earthquakes_spde102 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range10, .5), 
                                          prior.sigma = c(prior.median.sd102, .5))

s_index102 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde102$n.spde,
                                   n.group=12)

stack_est102 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index102) , tag="est")

stack_pred102 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred1),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid1)),
                                                    Time=GroupTime1, UTM_x=Cov1$UTM_x, 
                                                    UTM_y=Cov1$UTM_y, DEPTH=Cov1$DEPTH), s_index102), 
                            tag="pred")
stack102 <- inla.stack(stack_est102, stack_pred102) 

#===========6.5 Run INLA program
#===for 10km x 10km
#=========6.5.1 Range: 10km and Standard Deviation 0.1
formula1001 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde1001,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output1001 <- inla(formula1001,
                   data=inla.stack.data(stack1001, spde=Earthquakes_spde1001),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack1001), compute=TRUE))   
summary(output1001)

#=========6.5.2 Range: 10km and Standard Deviation 0.5
formula1005 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde1005,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output1005 <- inla(formula1005,
                   data=inla.stack.data(stack1005, spde=Earthquakes_spde1005),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack1005), compute=TRUE)) 
summary(output1005)

#=========6.5.3 Range: 10km and Standard Deviation 1
formula101 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde101,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output101 <- inla(formula101,
                  data=inla.stack.data(stack101, spde=Earthquakes_spde101),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack101), compute=TRUE))   
summary(output101)

#=========6.5.4 Range: 10km and Standard Deviation 1.5
formula1015 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde1015,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output1015 <- inla(formula1015,
                   data=inla.stack.data(stack1015, spde=Earthquakes_spde1015),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack1015), compute=TRUE)) 
summary(output1015)

#=========6.5.5 Range: 10km and Standard Deviation 2
formula102 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde102,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output102 <- inla(formula102,
                  data=inla.stack.data(stack102, spde=Earthquakes_spde102),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack102), compute=TRUE))   
summary(output102)

#======= Do the same steps like 10km x 10km for 20km, 30km, 40km, and 50km with following syntax
#=========================== 6.2 20km * 20km
Grid.CoordUTM2<-as.data.frame(coordsGrid.UTM2)
m2<-nrow(Grid.CoordUTM2)
Grid.CoordUTM2<-as.matrix(Grid.CoordUTM2)
GroupTime2<-rep(c(1:12),each=m2)
All.Grid2<-kronecker(matrix(1, 12, 1), Grid.CoordUTM2)

#============ 6.2.1 PENALIZED COMPLEXITY : 20km, 0.1
prior.median.sd2001 = 0.1; prior.median.range20 = 20000
Earthquakes_spde2001 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range20, .5), 
                                           prior.sigma = c(prior.median.sd2001, .5))
s_index2001 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde2001$n.spde,
                                    n.group=12)
stack_est2001 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index2001) , tag="est")
A_pred2 <- inla.spde.make.A(mesh=Earthquakes_mesh,
                            loc=as.matrix(All.Grid2),
                            group=GroupTime2,  #selected day for prediction
                            n.group=12)
Cov2<-data.frame(UTM_x=All.Grid2[,1]/1000,UTM_y=All.Grid2[,2]/1000,DEPTH=NA,Time=GroupTime2)
stack_pred2001 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred2),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid2)),
                                                     Time=GroupTime2, UTM_x=Cov2$UTM_x, 
                                                     UTM_y=Cov2$UTM_y, DEPTH=Cov2$DEPTH), s_index2001), 
                             tag="pred")
stack2001 <- inla.stack(stack_est2001, stack_pred2001) 
formula2001 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde2001,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output2001 <- inla(formula2001,
                   data=inla.stack.data(stack2001, spde=Earthquakes_spde2001),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack2001), compute=TRUE))   
summary(output2001)

#============ 7.6.2 PENALIZED COMPLEXITY : 20km, 0.5
prior.median.sd2005 = 0.5; prior.median.range20 = 20000
Earthquakes_spde2005 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range20, .5), 
                                           prior.sigma = c(prior.median.sd2005, .5))
s_index2005 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde2005$n.spde,
                                    n.group=12)
stack_est2005 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index2005) , tag="est")
stack_pred2005 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred2),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid2)),
                                                     Time=GroupTime2, UTM_x=Cov2$UTM_x, 
                                                     UTM_y=Cov2$UTM_y, DEPTH=Cov2$DEPTH), s_index2005), 
                             tag="pred")
stack2005 <- inla.stack(stack_est2005, stack_pred2005) 
formula2005 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde2005,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output2005 <- inla(formula2005,
                   data=inla.stack.data(stack2005, spde=Earthquakes_spde2005),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack2005), compute=TRUE))   
summary(output2005)

#============ 6.2.3 PENALIZED COMPLEXITY : 20km, 1
prior.median.sd201 = 1; prior.median.range20 = 20000
Earthquakes_spde201 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range20, .5), 
                                          prior.sigma = c(prior.median.sd201, .5))
s_index201 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde201$n.spde,
                                   n.group=12)
stack_est201 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index201) , tag="est")
stack_pred201 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred2),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid2)),
                                                    Time=GroupTime2, UTM_x=Cov2$UTM_x, 
                                                    UTM_y=Cov2$UTM_y, DEPTH=Cov2$DEPTH), s_index201), 
                            tag="pred")
stack201 <- inla.stack(stack_est201, stack_pred201) 
formula201 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde201,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output201 <- inla(formula201,
                  data=inla.stack.data(stack201, spde=Earthquakes_spde201),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack201), compute=TRUE))   
summary(output201)

#============ 6.2.4 PENALIZED COMPLEXITY : 20km, 1.5
prior.median.sd2015 = 1.5; prior.median.range20 = 20000
Earthquakes_spde2015 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range20, .5), 
                                           prior.sigma = c(prior.median.sd2015, .5))
s_index2015 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde2015$n.spde,
                                    n.group=12)
stack_est2015 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index2015) , tag="est")
stack_pred2015 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred2),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid2)),
                                                     Time=GroupTime2, UTM_x=Cov2$UTM_x, 
                                                     UTM_y=Cov2$UTM_y, DEPTH=Cov2$DEPTH), s_index2015), 
                             tag="pred")
stack2015 <- inla.stack(stack_est2015, stack_pred2015) 
formula2015 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde2015,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output2015 <- inla(formula2015,
                   data=inla.stack.data(stack2015, spde=Earthquakes_spde2015),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack2015), compute=TRUE))   
summary(output2015)

#============ 6.2.5 PENALIZED COMPLEXITY : 20km, 2
prior.median.sd202 = 2; prior.median.range20 = 20000
Earthquakes_spde202 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range20, .5), 
                                          prior.sigma = c(prior.median.sd202, .5))
s_index202 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde202$n.spde,
                                   n.group=12)
stack_est202 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index202) , tag="est")
stack_pred202 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred2),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid2)),
                                                    Time=GroupTime2, UTM_x=Cov2$UTM_x, 
                                                    UTM_y=Cov2$UTM_y, DEPTH=Cov2$DEPTH), s_index202), 
                            tag="pred")
stack202 <- inla.stack(stack_est202, stack_pred202) 
formula202 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde202,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output202 <- inla(formula202,
                  data=inla.stack.data(stack202, spde=Earthquakes_spde202),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack202), compute=TRUE))   
summary(output202)

#=========================== 6.3 30km * 30km
Grid.CoordUTM3<-as.data.frame(coordsGrid.UTM3)
m3<-nrow(Grid.CoordUTM3)
Grid.CoordUTM3<-as.matrix(Grid.CoordUTM3)
GroupTime3<-rep(c(1:12),each=m3)
All.Grid3<-kronecker(matrix(1, 12, 1), Grid.CoordUTM3)

#============ 6.3.1 PENALIZED COMPLEXITY : 30km, 0.1
prior.median.sd3001 = 0.1; prior.median.range30 = 30000
Earthquakes_spde3001 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range30, .5), 
                                           prior.sigma = c(prior.median.sd3001, .5))
s_index3001 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde3001$n.spde,
                                    n.group=12)
stack_est3001 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index3001) , tag="est")
A_pred3 <- inla.spde.make.A(mesh=Earthquakes_mesh,
                            loc=as.matrix(All.Grid3),
                            group=GroupTime3,  #selected day for prediction
                            n.group=12)
Cov3<-data.frame(UTM_x=All.Grid3[,1]/1000,UTM_y=All.Grid3[,2]/1000,DEPTH=NA,Time=GroupTime3)
stack_pred3001 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred3),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid3)),
                                                     Time=GroupTime3, UTM_x=Cov3$UTM_x, 
                                                     UTM_y=Cov3$UTM_y, DEPTH=Cov3$DEPTH), s_index3001), 
                             tag="pred")
stack3001 <- inla.stack(stack_est3001, stack_pred3001) 
formula3001 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde3001,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output3001 <- inla(formula3001,
                   data=inla.stack.data(stack3001, spde=Earthquakes_spde3001),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack3001), compute=TRUE))   
summary(output3001)

#============ 6.3.2 PENALIZED COMPLEXITY : 30km, 0.5
prior.median.sd3005 = 0.5; prior.median.range30 = 30000
Earthquakes_spde3005 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range30, .5), 
                                           prior.sigma = c(prior.median.sd3005, .5))
s_index3005 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde3005$n.spde,
                                    n.group=12)
stack_est3005 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index3005) , tag="est")
stack_pred3005 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred3),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid3)),
                                                     Time=GroupTime3, UTM_x=Cov3$UTM_x, 
                                                     UTM_y=Cov3$UTM_y, DEPTH=Cov3$DEPTH), s_index3005), 
                             tag="pred")
stack3005 <- inla.stack(stack_est3005, stack_pred3005) 
formula3005 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde3005,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output3005 <- inla(formula3005,
                   data=inla.stack.data(stack3005, spde=Earthquakes_spde3005),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack3005), compute=TRUE))   
summary(output3005)

#============ 6.3.3 PENALIZED COMPLEXITY : 30km, 1
prior.median.sd301 = 1; prior.median.range30 = 30000
Earthquakes_spde301 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range30, .5), 
                                          prior.sigma = c(prior.median.sd301, .5))
s_index301 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde301$n.spde,
                                   n.group=12)
stack_est301 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index301) , tag="est")
stack_pred301 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred3),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid3)),
                                                    Time=GroupTime3, UTM_x=Cov3$UTM_x, 
                                                    UTM_y=Cov3$UTM_y, DEPTH=Cov3$DEPTH), s_index301), 
                            tag="pred")
stack301 <- inla.stack(stack_est301, stack_pred301) 
formula301 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde301,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output301 <- inla(formula301,
                  data=inla.stack.data(stack301, spde=Earthquakes_spde301),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack301), compute=TRUE))   
summary(output301)

#============ 6.3.4 PENALIZED COMPLEXITY : 30km, 1.5
prior.median.sd3015 = 1.5; prior.median.range30 = 30000
Earthquakes_spde3015 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range30, .5), 
                                           prior.sigma = c(prior.median.sd3015, .5))
s_index3015 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde3015$n.spde,
                                    n.group=12)
stack_est3015 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index3015) , tag="est")
stack_pred3015 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred3),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid3)),
                                                     Time=GroupTime3, UTM_x=Cov3$UTM_x, 
                                                     UTM_y=Cov3$UTM_y, DEPTH=Cov3$DEPTH), s_index3015), 
                             tag="pred")
stack3015 <- inla.stack(stack_est3015, stack_pred3015) 
formula3015 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde3015,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output3015 <- inla(formula3015,
                   data=inla.stack.data(stack3015, spde=Earthquakes_spde3015),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack3015), compute=TRUE))   
summary(output3015)

#============ 6.3.5 PENALIZED COMPLEXITY : 30km, 2
prior.median.sd302 = 2; prior.median.range30 = 30000
Earthquakes_spde302 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range30, .5), 
                                          prior.sigma = c(prior.median.sd302, .5))
s_index302 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde302$n.spde,
                                   n.group=12)
stack_est302 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index302) , tag="est")
stack_pred302 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred3),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid3)),
                                                    Time=GroupTime3, UTM_x=Cov3$UTM_x, 
                                                    UTM_y=Cov3$UTM_y, DEPTH=Cov3$DEPTH), s_index302), 
                            tag="pred")
stack302 <- inla.stack(stack_est302, stack_pred302) 
formula302 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde302,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output302 <- inla(formula302,
                  data=inla.stack.data(stack302, spde=Earthquakes_spde302),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack302), compute=TRUE))   
summary(output302)

#=========================== 6.4 40km * 40km
Grid.CoordUTM4<-as.data.frame(coordsGrid.UTM4)
m4<-nrow(Grid.CoordUTM4)
Grid.CoordUTM4<-as.matrix(Grid.CoordUTM4)
GroupTime4<-rep(c(1:12),each=m4)
All.Grid4<-kronecker(matrix(1, 12, 1), Grid.CoordUTM4)

#============ 6.4.1 PENALIZED COMPLEXITY : 40km, 0.1
prior.median.sd4001 = 0.1; prior.median.range40 = 40000
Earthquakes_spde4001 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range40, .5), 
                                           prior.sigma = c(prior.median.sd4001, .5))
s_index4001 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde4001$n.spde,
                                    n.group=12)
stack_est4001 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index4001) , tag="est")
A_pred4 <- inla.spde.make.A(mesh=Earthquakes_mesh,
                            loc=as.matrix(All.Grid4),
                            group=GroupTime4,  #selected day for prediction
                            n.group=12)
Cov4<-data.frame(UTM_x=All.Grid4[,1]/1000,UTM_y=All.Grid4[,2]/1000,DEPTH=NA,Time=GroupTime4)
stack_pred4001 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred4),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid4)),
                                                     Time=GroupTime4, UTM_x=Cov4$UTM_x, 
                                                     UTM_y=Cov4$UTM_y, DEPTH=Cov4$DEPTH), s_index4001), 
                             tag="pred")
stack4001 <- inla.stack(stack_est4001, stack_pred4001) 
formula4001 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde4001,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output4001 <- inla(formula4001,
                   data=inla.stack.data(stack4001, spde=Earthquakes_spde4001),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack4001), compute=TRUE))   
summary(output4001)

#============ 6.4.2 PENALIZED COMPLEXITY : 40km, 0.5
prior.median.sd4005 = 0.5; prior.median.range40 = 40000
Earthquakes_spde4005 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range40, .5), 
                                           prior.sigma = c(prior.median.sd4005, .5))
s_index4005 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde4005$n.spde,
                                    n.group=12)
stack_est4005 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index4005) , tag="est")
stack_pred4005 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred4),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid4)),
                                                     Time=GroupTime4, UTM_x=Cov4$UTM_x, 
                                                     UTM_y=Cov4$UTM_y, DEPTH=Cov4$DEPTH), s_index4005), 
                             tag="pred")
stack4005 <- inla.stack(stack_est4005, stack_pred4005) 
formula4005 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde4005,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output4005 <- inla(formula4005,
                   data=inla.stack.data(stack4005, spde=Earthquakes_spde4005),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack4005), compute=TRUE))   
summary(output4005)

#============ 6.4.3 PENALIZED COMPLEXITY : 40km, 1
prior.median.sd401 = 1; prior.median.range40 = 40000
Earthquakes_spde401 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range40, .5), 
                                          prior.sigma = c(prior.median.sd401, .5))
s_index401 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde401$n.spde,
                                   n.group=12)
stack_est401 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index401) , tag="est")
stack_pred401 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred4),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid4)),
                                                    Time=GroupTime4, UTM_x=Cov4$UTM_x, 
                                                    UTM_y=Cov4$UTM_y, DEPTH=Cov4$DEPTH), s_index401), 
                            tag="pred")
stack401 <- inla.stack(stack_est401, stack_pred401) 
formula401 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde401,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output401 <- inla(formula401,
                  data=inla.stack.data(stack401, spde=Earthquakes_spde401),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack401), compute=TRUE))   
summary(output401)

#============ 6.4.4 PENALIZED COMPLEXITY : 40km, 1.5
prior.median.sd4015 = 1.5; prior.median.range40 = 40000
Earthquakes_spde4015 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range40, .5), 
                                           prior.sigma = c(prior.median.sd4015, .5))
s_index4015 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde4015$n.spde,
                                    n.group=12)
stack_est4015 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index4015) , tag="est")
stack_pred4015 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred4),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid4)),
                                                     Time=GroupTime4, UTM_x=Cov4$UTM_x, 
                                                     UTM_y=Cov4$UTM_y, DEPTH=Cov4$DEPTH), s_index4015), 
                             tag="pred")
stack4015 <- inla.stack(stack_est4015, stack_pred4015) 
formula4015 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde4015,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output4015 <- inla(formula4015,
                   data=inla.stack.data(stack4015, spde=Earthquakes_spde4015),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack4015), compute=TRUE))   
summary(output4015)

#============ 6.4.5 PENALIZED COMPLEXITY : 40km, 2
prior.median.sd402 = 2; prior.median.range40 = 40000
Earthquakes_spde402 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range40, .5), 
                                          prior.sigma = c(prior.median.sd402, .5))
s_index402 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde402$n.spde,
                                   n.group=12)
stack_est402 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index402) , tag="est")
stack_pred402 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred4),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid4)),
                                                    Time=GroupTime4, UTM_x=Cov4$UTM_x, 
                                                    UTM_y=Cov4$UTM_y, DEPTH=Cov4$DEPTH), s_index402), 
                            tag="pred")
stack402 <- inla.stack(stack_est402, stack_pred402) 
formula402 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde402,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output402 <- inla(formula402,
                  data=inla.stack.data(stack402, spde=Earthquakes_spde402),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack402), compute=TRUE))   
summary(output402)

#=========================== 6.5 50km * 50km
Grid.CoordUTM5<-as.data.frame(coordsGrid.UTM5)
m5<-nrow(Grid.CoordUTM5)
Grid.CoordUTM5<-as.matrix(Grid.CoordUTM5)
GroupTime5<-rep(c(1:12),each=m5)
All.Grid5<-kronecker(matrix(1, 12, 1), Grid.CoordUTM5)

#============ 6.5.1 PENALIZED COMPLEXITY : 50km, 0.1
prior.median.sd5001 = 0.1; prior.median.range50 = 50000
Earthquakes_spde5001 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range50, .5), 
                                           prior.sigma = c(prior.median.sd5001, .5))
s_index5001 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde5001$n.spde,
                                    n.group=12)
stack_est5001 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index5001) , tag="est")
A_pred5 <- inla.spde.make.A(mesh=Earthquakes_mesh,
                            loc=as.matrix(All.Grid5),
                            group=GroupTime5,  #selected day for prediction
                            n.group=12)
Cov5<-data.frame(UTM_x=All.Grid5[,1]/1000,UTM_y=All.Grid5[,2]/1000,DEPTH=NA,Time=GroupTime5)
stack_pred5001 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred5),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid5)),
                                                     Time=GroupTime5, UTM_x=Cov5$UTM_x, 
                                                     UTM_y=Cov5$UTM_y, DEPTH=Cov5$DEPTH), s_index5001), 
                             tag="pred")
stack5001 <- inla.stack(stack_est5001, stack_pred5001) 
formula5001 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde5001,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output5001 <- inla(formula5001,
                   data=inla.stack.data(stack5001, spde=Earthquakes_spde5001),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack5001), compute=TRUE))   
summary(output5001)

#============ 6.5.2 PENALIZED COMPLEXITY : 50km, 0.5
prior.median.sd5005 = 0.5; prior.median.range50 = 50000
Earthquakes_spde5005 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range50, .5), 
                                           prior.sigma = c(prior.median.sd5005, .5))
s_index5005 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde5005$n.spde,
                                    n.group=12)
stack_est5005 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index5005) , tag="est")
stack_pred5005 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred5),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid5)),
                                                     Time=GroupTime5, UTM_x=Cov5$UTM_x, 
                                                     UTM_y=Cov5$UTM_y, DEPTH=Cov5$DEPTH), s_index5005), 
                             tag="pred")
stack5005 <- inla.stack(stack_est5005, stack_pred5005) 
formula5005 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde5005,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output5005 <- inla(formula5005,
                   data=inla.stack.data(stack5005, spde=Earthquakes_spde5005),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack5005), compute=TRUE))   
summary(output5005)

#============ 6.5.3 PENALIZED COMPLEXITY : 50km, 1
prior.median.sd501 = 1; prior.median.range50 = 50000
Earthquakes_spde501 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range50, .5), 
                                          prior.sigma = c(prior.median.sd501, .5))
s_index501 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde501$n.spde,
                                   n.group=12)
stack_est501 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index501) , tag="est")
stack_pred501 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred5),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid5)),
                                                    Time=GroupTime5, UTM_x=Cov5$UTM_x, 
                                                    UTM_y=Cov5$UTM_y, DEPTH=Cov5$DEPTH), s_index501), 
                            tag="pred")
stack501 <- inla.stack(stack_est501, stack_pred501) 
formula501 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde501,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output501 <- inla(formula501,
                  data=inla.stack.data(stack501, spde=Earthquakes_spde501),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack501), compute=TRUE))   
summary(output501)

#============ 6.5.4 PENALIZED COMPLEXITY : 50km, 1.5
prior.median.sd5015 = 1.5; prior.median.range50 = 50000
Earthquakes_spde5015 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                           prior.range = c(prior.median.range50, .5), 
                                           prior.sigma = c(prior.median.sd5015, .5))
s_index5015 <- inla.spde.make.index(name="spatial.field",
                                    n.spde=Earthquakes_spde5015$n.spde,
                                    n.group=12)
stack_est5015 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                            A=list(1,A_est),
                            effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                    Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                    UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                         s_index5015) , tag="est")
stack_pred5015 <- inla.stack(data=list(MAG=NA),
                             A=list(1, A_pred5),
                             effects=list(data.frame(Intercept = rep(1, nrow(All.Grid5)),
                                                     Time=GroupTime5, UTM_x=Cov5$UTM_x, 
                                                     UTM_y=Cov5$UTM_y, DEPTH=Cov5$DEPTH), s_index5015), 
                             tag="pred")
stack5015 <- inla.stack(stack_est5015, stack_pred5015) 

formula5015 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde5015,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output5015 <- inla(formula5015,
                   data=inla.stack.data(stack5015, spde=Earthquakes_spde5015),
                   family="gaussian",control.compute = control$compute,
                   control.predictor=list(A=inla.stack.A(stack5015), compute=TRUE))   
summary(output5015)

#============ 6.5.5 PENALIZED COMPLEXITY : 50km, 2
prior.median.sd502 = 2; prior.median.range50 = 50000
Earthquakes_spde502 = inla.spde2.pcmatern(mesh=Earthquakes_mesh, 
                                          prior.range = c(prior.median.range50, .5), 
                                          prior.sigma = c(prior.median.sd502, .5))
s_index502 <- inla.spde.make.index(name="spatial.field",
                                   n.spde=Earthquakes_spde502$n.spde,
                                   n.group=12)
stack_est502 <- inla.stack(data=list(MAG=log(Earthquakes20$MAG)),
                           A=list(1,A_est),
                           effects=list(data.frame(Intercept = rep(1, nrow(Earthquakes20)), 
                                                   Time=Earthquakes20$Month, UTM_x=Earthquakes20$UTM_x/1000, 
                                                   UTM_y=Earthquakes20$UTM_y/1000, DEPTH=Earthquakes20$DEPTH), 
                                        s_index502) , tag="est")
stack_pred502 <- inla.stack(data=list(MAG=NA),
                            A=list(1, A_pred5),
                            effects=list(data.frame(Intercept = rep(1, nrow(All.Grid5)),
                                                    Time=GroupTime5, UTM_x=Cov5$UTM_x, 
                                                    UTM_y=Cov5$UTM_y, DEPTH=Cov5$DEPTH), s_index502), 
                            tag="pred")
stack502 <- inla.stack(stack_est502, stack_pred502) 
formula502 <- MAG ~ -1 + Intercept + UTM_x + UTM_y + 
  f(spatial.field, model=Earthquakes_spde502,group=spatial.field.group, 
    control.group=list(model="rw1"))

# ATTENTION: the run is computationally intensive!
output502 <- inla(formula502,
                  data=inla.stack.data(stack502, spde=Earthquakes_spde502),
                  family="gaussian",control.compute = control$compute,
                  control.predictor=list(A=inla.stack.A(stack502), compute=TRUE))   
summary(output502)

#========================================================================================================================
#==================== 7. Model Selection
#========================================================================================================================
#===========7.1 DIC
range_km <- c(10,20,30,40,50,
              10,20,30,40,50,
              10,20,30,40,50,
              10,20,30,40,50,
              10,20,30,40,50)
std <- c(0.1, 0.5, 1, 1.5, 2,
         0.1, 0.5, 1, 1.5, 2,
         0.1, 0.5, 1, 1.5, 2,
         0.1, 0.5, 1, 1.5, 2,
         0.1, 0.5, 1, 1.5, 2)
DIC2 <- c(7.84, 24.57, 31.00, 32.27, 33.24, 
          8.30, 27.07, 30.96, 31.19, 32.78, 
          7.71, 25.93, 30.83, 32.37, 33.00, 
          4.25, 25.99, 30.87, 31.95, 32.52, 
          6.21, 29.66, 30.77, 31.45, 32.73)
z2 <- matrix(DIC2, nrow = length(unique(std)), byrow = TRUE)

plot_ly(x = ~std, y = ~range_km, z = ~z2, type = "surface", 
        colorscale = list(list(0, "cornsilk"), list(0.5, "bisque"), list(1, "darkorange")),
        colorbar = list(title = "DIC", orientation = "h", x = 0.5, y = -0.2)) %>%
  layout(title = "Surface Plot of DIC",
         scene = list(xaxis = list(title = "Standard Deviation"),
                      yaxis = list(title = "Range (km)"),
                      zaxis = list(title = "DIC")))

#===========7.2 WAIC
WAIC2 <- c(107.91,151.60,146.56,150.09,149.62,
           114.72,140.67,146.32,153.23,150.63,
           106.04,144.30,146.92,148.32,149.80,
           110.05,149.58,145.92,149.65,151.11,
           109.32,133.52,145.93,150.44,150.41)
z2 <- matrix(WAIC2, nrow = length(unique(std)), byrow = TRUE)

plot_ly(x = ~std, y = ~range_km, z = ~z2, type = "surface", 
        colorscale = list(list(0, "cornsilk"), list(0.5, "bisque"), list(1, "darkorange2")),
        colorbar = list(title = "DIC", orientation = "h", x = 0.5, y = -0.2)) %>%
  layout(title = "Surface Plot of DIC",
         scene = list(xaxis = list(title = "Standard Deviation"),
                      yaxis = list(title = "Range (km)"),
                      zaxis = list(title = "WAIC")))

#========================================================================================================================
#==================== 8. Exceedance Probability
#========================================================================================================================
MAG<-log(5)
ExProb <- unlist(lapply(output1001$marginals.fitted.values, function(X) {
  if (all(is.finite(X))) {
    return(1 - inla.pmarginal(MAG, X))
  } else {
    return(NA)
  }
}))

#========================================================================================================================
#==================== 9. Visualization
#========================================================================================================================
#=========== 9.1 Prepare output visualization in map
index_pred <- inla.stack.index(stack1001,"pred")$data
index_est <- inla.stack.index(stack1001,"est")$data

Earthquakes_Pred<-data.frame(x=All.Grid1[,1],y=All.Grid1[,2],Time=GroupTime1) 
Earthquakes_Est<-data.frame(x=Earthquakes20$UTM_x,y=Earthquakes20$UTM_y, Time=Earthquakes20$Month) 

Earthquakes_Pred$pred_mean <- exp(output1001$summary.fitted.values[index_pred, "mean"])
Earthquakes_Pred$pred_ll <- exp(output1001$summary.fitted.values[index_pred, "0.025quant"])
Earthquakes_Pred$pred_ul <- exp(output1001$summary.fitted.values[index_pred, "0.975quant"])

Earthquakes_Est$pred_mean <- exp(output1001$summary.fitted.values[index_est, "mean"])
Earthquakes_Est$pred_ll <- exp(output1001$summary.fitted.values[index_est, "0.025quant"])
Earthquakes_Est$pred_ul <- exp(output1001$summary.fitted.values[index_est, "0.975quant"])
summary(Earthquakes_Est$pred_mean)

Earthquakes_Values<-rbind(Earthquakes_Est,Earthquakes_Pred)
Earthquakes_Pred$pred_prob <- ExProb[index_pred]

write.csv(Earthquakes_Pred,"Earthquakes_Pred.csv")
write.csv(Earthquakes_Est,"Earthquakes_Est.csv")
Earthquakes_Pred<-read.csv("Earthquakes_Pred.csv", sep=",")

summary(Earthquakes_Est$pred_mean)
summary(Earthquakes_Est$pred_mean)

#=========== 9.2 High-Resolution of Earthquake Interpolation
#========== 9.2.1 Overall: By Month
Earthquakes_Pred$title <- "By Month"
all = ggplot() + 
  geom_tile(data = Earthquakes_Pred, aes(x = x, y = y, fill = pred_mean)) + 
  geom_contour(data = Earthquakes_Pred, aes(x = x, y = y, z = pred_mean, colour = stat(level)), size = 0.2, colour = "red", alpha = 0.5) +
  scale_fill_gradientn(colours = c("gray99", "yellow", "red")) +
  facet_wrap(~ Time, labeller = as_labeller(Month_names), ncol = 4) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  geom_sf(data = Sumatra_sf_utm1, fill = NA, size=0.01) +  
  theme(legend.position = "bottom", text = element_text(size = 19)) + 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1,  
                               frame.colour = "black",
                               barwidth = 20,
                               barheight = 1.5)) +
  theme(panel.background = element_rect(fill = "white", color = NA), 
        panel.grid.major = element_line(color = "grey95", linetype = "solid"),
        panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
        axis.text.x = element_blank(), #remove x axis labels
        axis.ticks.x = element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),  #remove y axis labels
        axis.ticks.y = element_blank()  #remove y axis ticks
  ) + 
  labs(fill = "Earthquakes \nMagnitude")
all

#========== 9.2.2 Overall: By Province
Earthquakes_Pred_dec <- Earthquakes_Pred[Earthquakes_Pred$Time == 12, ]
Earthquakes_Pred_dec$title = "By Province"
ggplot() +    
  geom_tile(data = Earthquakes_Pred_dec, aes(x = x, y = y, fill = pred_mean)) +    
  geom_contour(data = Earthquakes_Pred_dec, aes(x = x, y = y, z = pred_mean, colour = stat(level)), 
               size = 0.2, colour = "red", alpha = 0.5) +   
  scale_fill_gradientn(colours = c("gray99", "yellow", "red")) +   
  theme_bw() +    
  ylab("Latitude") +    
  xlab("Longitude") +    
  geom_sf(data = Sumatra_sf_utm1, fill = NA, size = 0.01) +     
  geom_sf_text(data =  Sumatra_sf_utm1, aes(label = NAME_1), size = 3, color = "black", fontface = "bold") +
  geom_sf(data = patahan_sumatera, color = "blue", size = 1) +
  theme(legend.position = "bottom", 
        text = element_text(size = 12)) +   
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
    plot.title = element_text(hjust = 0.5, size = 12),   # Menyesuaikan ukuran dan posisi judul
    panel.grid = element_blank()                         # Menghilangkan grid
  )+
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1, 
                               frame.colour = "black", 
                               barwidth = 20, 
                               barheight = 1.5))+
  facet_grid(. ~ title)

#========== 9.2.3 District Level
max_pred = max(Earthquakes_Pred$pred_mean)
kab = st_read(file.choose())
kab <- st_transform(kab, crs = st_crs(Sumatra))
kab_valid <- st_make_valid(kab)

#Bengkulu Province
Sumatra_Bengkulu <- Sumatra_sf_utm1 %>% filter(NAME_1 == 'Bengkulu')
Sumatra_Bengkulu_wgs84 <- st_transform(Sumatra_Bengkulu, crs = st_crs(kab))
bbox_Sumatra_Bengkulu <- st_bbox(Sumatra_Bengkulu_wgs84)
kab_cropped <- st_crop(kab_valid, bbox_Sumatra_Bengkulu)
kab_bengkulu<-kab_cropped[kab_cropped$KAB_KOTA %in%c("MUKOMUKO","LEBONG","BENGKULU UTARA", "REJANG LEBONG", "BENGKULU TENGAH", "KEPAHIANG","KOTA BENGKULU", "SELUMA", "BENGKULU SELATAN", "KAUR"),]
kab_bengkulu_sf <- st_as_sf(kab_bengkulu)
kab_bengkulu_sf_utm1 <- st_transform(kab_bengkulu_sf, crs = 32648)

Earthquakes_Pred_Bengkulu <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kab_bengkulu_sf_utm1)$xmin) & x <= max(st_bbox(kab_bengkulu_sf_utm1)$xmax) &
           y >= min(st_bbox(kab_bengkulu_sf_utm1)$ymin) & y <= max(st_bbox(kab_bengkulu_sf_utm1)$ymax))

min(Earthquakes_Pred_Bengkulu$pred_mean)
max(Earthquakes_Pred_Bengkulu$pred_mean)

#==== plot with smoothness
Earthquakes_Pred_Bengkulu$x <- jitter(Earthquakes_Pred_Bengkulu$x, amount = 0.0001)
Earthquakes_Pred_Bengkulu$y <- jitter(Earthquakes_Pred_Bengkulu$y, amount = 0.0001)

# Interpolation after jitter
interpolated_bengkulu <- with(Earthquakes_Pred_Bengkulu, 
                              interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

# Dataframe format
interpolated_df_beng <- as.data.frame(expand.grid(x = interpolated_bengkulu$x, y = interpolated_bengkulu$y))
interpolated_df_beng$pred_mean <- as.vector(interpolated_bengkulu$z)
interpolated_df_beng <- na.omit(interpolated_df_beng)  # Menghapus nilai NA

#Plot
interpolated_df_beng$title = "Bengkulu"
bengkulu = ggplot() + 
  geom_raster(data = interpolated_df_beng, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kab_bengkulu_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kab_bengkulu, aes(label = KAB_KOTA), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
bengkulu

#West Sumatra Province
Sumatra_Barat <- Sumatra_sf_utm1 %>% filter(NAME_1 == 'Sumatera Barat')
Sumatra_Barat_wgs84 <- st_transform(Sumatra_Barat, crs = st_crs(kab))
bbox_Sumatra_Barat <- st_bbox(Sumatra_Barat_wgs84)
kab_cropped <- st_crop(kab_valid, bbox_Sumatra_Barat)
kab_sumbar<-kab_cropped[kab_cropped$KAB_KOTA %in%c("AGAM","DHARMASRAYA", "KEPULAUAN MENTAWAI", "LIMA PULUH KOTA","PADANG PARIAMAN", "PASAMAN", "PASAMAN BARAT", "PESISIR SELATAN","SIJUNJUNG", "SOLOK", "SOLOK SELATAN", "TANAH DATAR", "KOTA BUKITTINGGI", "KOTA PADANG", "KOTA PADANG PANJANG", "KOTA PARIAMAN","KOTA PAYAKUMBUH",  "KOTA SAWAHLUNTO", "KOTA SOLOK"),]
kab_sumbar_sf <- st_as_sf(kab_sumbar)
kab_sumbar_sf_utm1 <- st_transform(kab_sumbar_sf, crs = 32648)

Earthquakes_Pred_Sumbar <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kab_sumbar_sf_utm1)$xmin) & x <= max(st_bbox(kab_sumbar_sf_utm1)$xmax) &
           y >= min(st_bbox(kab_sumbar_sf_utm1)$ymin) & y <= max(st_bbox(kab_sumbar_sf_utm1)$ymax))
min(Earthquakes_Pred_Sumbar$pred_mean)
max(Earthquakes_Pred_Sumbar$pred_mean)

Earthquakes_Pred_Sumbar$x <- jitter(Earthquakes_Pred_Sumbar$x, amount = 0.0001)
Earthquakes_Pred_Sumbar$y <- jitter(Earthquakes_Pred_Sumbar$y, amount = 0.0001)

interpolated_sumbar <- with(Earthquakes_Pred_Sumbar, 
                            interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_sumbar <- as.data.frame(expand.grid(x = interpolated_sumbar$x, y = interpolated_sumbar$y))
interpolated_df_sumbar$pred_mean <- as.vector(interpolated_sumbar$z)
interpolated_df_sumbar <- na.omit(interpolated_df_sumbar)  # Menghapus nilai NA

#Plot
interpolated_df_sumbar$title = "Sumatra Barat"
sumbar = ggplot() + 
  geom_raster(data = interpolated_df_sumbar, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kab_sumbar_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kab_sumbar, aes(label = KAB_KOTA), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
sumbar

#Nort Sumatra Province
Sumatra_Utara <- Sumatra_sf_utm1 %>% filter(NAME_1 == 'Sumatera Utara')
Sumatra_Utara_wgs84 <- st_transform(Sumatra_Utara, crs = st_crs(kab))
bbox_Sumatra_Utara <- st_bbox(Sumatra_Utara_wgs84)
kab_cropped <- st_crop(kab_valid, bbox_Sumatra_Utara)

kab_sumut<-kab_cropped[kab_cropped$KAB_KOTA %in%c("ASAHAN", "BATU BARA", "DAIRI", "DELI SERDANG", "HUMBANG HASUNDUTAN","KARO", "LABUHANBATU", "LABUHANBATU SELATAN", "LABUHANBATU UTARA", "LANGKAT","MANDAILING NATAL", "NIAS", "NIAS BARAT","NIAS SELATAN", "NIAS UTARA", "PADANG LAWAS","PADANG LAWAS UTARA","PAKPAK BHARAT", "SAMOSIR","SERDANG BEDAGAI", "SIMALUNGUN","TAPANULI SELATAN","TAPANULI TENGAH","TAPANULI UTARA","TOBA SAMOSIR","KOTA BINJAI","KOTA GUNUNGSITOLI", "KOTA MEDAN", "KOTA PADANG SIDIMPUAN","KOTA PEMATANGSIANTAR","KOTA SIBOLGA","KOTA TANJUNG BALAI",
                                                  "KOTA TEBING TINGGI"),]
kab_sumut_sf <- st_as_sf(kab_sumut)
kab_sumut_sf_utm1 <- st_transform(kab_sumut_sf, crs = 32648)

Earthquakes_Pred_Sumut <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kab_sumut_sf_utm1)$xmin) & x <= max(st_bbox(kab_sumut_sf_utm1)$xmax) &
           y >= min(st_bbox(kab_sumut_sf_utm1)$ymin) & y <= max(st_bbox(kab_sumut_sf_utm1)$ymax))
#min(Earthquakes_Pred_Sumut$pred_mean)
#max(Earthquakes_Pred_Sumut$pred_mean)

Earthquakes_Pred_Sumut$x <- jitter(Earthquakes_Pred_Sumut$x, amount = 0.0001)
Earthquakes_Pred_Sumut$y <- jitter(Earthquakes_Pred_Sumut$y, amount = 0.0001)

interpolated_sumut <- with(Earthquakes_Pred_Sumut, 
                           interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_sumut <- as.data.frame(expand.grid(x = interpolated_sumut$x, y = interpolated_sumut$y))
interpolated_df_sumut$pred_mean <- as.vector(interpolated_sumut$z)
interpolated_df_sumut <- na.omit(interpolated_df_sumut)  # Menghapus nilai NA

#Plot
interpolated_df_sumut$title = "Sumatra Utara"
sumut = ggplot() + 
  geom_raster(data = interpolated_df_sumut, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kab_sumut_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kab_sumut, aes(label = KAB_KOTA), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
sumut

#========== 9.2.4 Sub-District Level
kec = st_read(file.choose())
kec <- st_transform(kec, crs = st_crs(Sumatra))
kec_valid <- st_make_valid(kec)

#=== Mukomuko, bengkulu utara, and lebong
kab_bengkulu_wgs84 <- st_transform(kab_bengkulu, crs = st_crs(kec))
bbox_kab_bengkulu <- st_bbox(kab_bengkulu_wgs84)
kec_cropped1 <- st_crop(kec_valid, bbox_kab_bengkulu)

kec_mbul<-kec_cropped1[kec_cropped1$KECAMATAN %in%c("V KOTO", "XIV KOTO","AIR DIKIT","AIR MAJUNTO","AIR RAMI", "IPUH","KOTA MUKOMUKO","LUBUK PINANG","MALIN DEMAN","PENARIK","PONDOK SUGUH","SELAGAN RAYA","SUNGAI RUMBAI","TERAMANG JAYA", "TERAS TERUNJAM", "AIR BESI","AIR NAPAL", "AIR PADANG", "ARMA JAYA","BATIK NAU", "GIRI MULYA","HULU PALIK","KERKAP","KETAHUN",
                                                    "KOTA ARGA MAKMUR","LAIS","MARGA SAKTI SEBELAT","NAPAL PUTIH","PADANG JAYA","PINANG RAYA","PUTRI HIJAU","TANJUNG AGUNG PALIK","ULOK KUPAI", "AMEN","BINGIN KUNING", "LEBONG ATAS", "LEBONG SAKTI","LEBONG SELATAN","LEBONG TENGAH","LEBONG UTARA","PINANG BELAPIS","RIMBO PENGADANG","TOPOS","TUBEI","URAM JAYA"),]
kec_mbul_sf <- st_as_sf(kec_mbul)
kec_mbul_sf_utm1 <- st_transform(kec_mbul_sf, crs = 32648)

Earthquakes_Pred_MBUL <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_mbul_sf_utm1)$xmin) & x <= max(st_bbox(kec_mbul_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_mbul_sf_utm1)$ymin) & y <= max(st_bbox(kec_mbul_sf_utm1)$ymax))
min(Earthquakes_Pred_MBUL$pred_mean)
max(Earthquakes_Pred_MBUL$pred_mean)

Earthquakes_Pred_MBUL$x <- jitter(Earthquakes_Pred_MBUL$x, amount = 0.0001)
Earthquakes_Pred_MBUL$y <- jitter(Earthquakes_Pred_MBUL$y, amount = 0.0001)

interpolated_mbul <- with(Earthquakes_Pred_MBUL, 
                          interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_mbul <- as.data.frame(expand.grid(x = interpolated_mbul$x, y = interpolated_mbul$y))
interpolated_df_mbul$pred_mean <- as.vector(interpolated_mbul$z)
interpolated_df_mbul <- na.omit(interpolated_df_mbul)

#Plot
interpolated_df_mbul$title = "Mukomuko, Bengkulu Utara, & Lebong"
mbul = ggplot() + 
  geom_raster(data = interpolated_df_mbul, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_mbul_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_mbul, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
mbul 

#=== Bengkulu Selatan and Kaur
kec_bsk<-kec_cropped1[kec_cropped1$KECAMATAN %in%c("KAUR SELATAN","KAUR UTARA","KELAM TENGAH",
                                                   "KAUR TENGAH","KINAL","LUNGKANG KULE","LUAS",
                                                   "MAJE","MUARA SAHUNG","NASAL","PADANG GUCI HILIR",
                                                   "PADANG GUCI HULU","SEMIDANG GUMAY","TANJUNG KEMUNING","TETAP",
                                                   "KEDURANG","SEGINIM","PINO","MANNA","KOTA MANNA","PINO RAYA",
                                                   "KEDURANG ILIR", "AIR NIPIS","ULU MANNA","BUNGA MAS",
                                                   "PASAR MANNA"),]
kec_bsk_sf <- st_as_sf(kec_bsk)
kec_bsk_sf_utm1 <- st_transform(kec_bsk_sf, crs = 32648)

Earthquakes_Pred_BSK <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_bsk_sf_utm1)$xmin) & x <= max(st_bbox(kec_bsk_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_bsk_sf_utm1)$ymin) & y <= max(st_bbox(kec_bsk_sf_utm1)$ymax))
min(Earthquakes_Pred_BSK$pred_mean)
max(Earthquakes_Pred_BSK$pred_mean)

Earthquakes_Pred_BSK$x <- jitter(Earthquakes_Pred_BSK$x, amount = 0.0001)
Earthquakes_Pred_BSK$y <- jitter(Earthquakes_Pred_BSK$y, amount = 0.0001)

interpolated_BSK <- with(Earthquakes_Pred_BSK, 
                         interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_BSK <- as.data.frame(expand.grid(x = interpolated_BSK$x, y = interpolated_BSK$y))
interpolated_df_BSK$pred_mean <- as.vector(interpolated_BSK$z)
interpolated_df_BSK <- na.omit(interpolated_df_BSK)  # Menghapus nilai NA

#Plot
interpolated_df_BSK$title = "Bengkulu Selatan and Kaur"
bsk = ggplot() + 
  geom_raster(data = interpolated_df_BSK, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_bsk_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_bsk, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
bsk

#=== Kab Mentawai
kab_sumbar_wgs84 <- st_transform(kab_sumbar, crs = st_crs(kec))
bbox_kab_sumbar <- st_bbox(kab_sumbar_wgs84)
kec_cropped <- st_crop(kec_valid, bbox_kab_sumbar)
kec_mentawai<-kec_cropped[kec_cropped$KECAMATAN %in%c("PAGAI SELATAN","PAGAI UTARA","SIBERUT BARAT","SIBERUT BARAT DAYA","SIBERUT SELATAN", "SIBERUT UTARA","SIBERUT TENGAH","SIKAKAP","SIPORA SELATAN","SIPORA UTARA"),]
kec_mentawai_sf <- st_as_sf(kec_mentawai)
kec_mentawai_sf_utm1 <- st_transform(kec_mentawai_sf, crs = 32648)

Earthquakes_Pred_Mentawai <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_mentawai_sf_utm1)$xmin) & x <= max(st_bbox(kec_mentawai_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_mentawai_sf_utm1)$ymin) & y <= max(st_bbox(kec_mentawai_sf_utm1)$ymax))

Earthquakes_Pred_Mentawai$x <- jitter(Earthquakes_Pred_Mentawai$x, amount = 0.0001)
Earthquakes_Pred_Mentawai$y <- jitter(Earthquakes_Pred_Mentawai$y, amount = 0.0001)

interpolated_mentawai <- with(Earthquakes_Pred_Mentawai, 
                              interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_mentawai <- as.data.frame(expand.grid(x = interpolated_mentawai$x, y = interpolated_mentawai$y))
interpolated_df_mentawai$pred_mean <- as.vector(interpolated_mentawai$z)
interpolated_df_mentawai <- na.omit(interpolated_df_mentawai)  

interpolated_df_mentawai$title = "Kep. Mentawai"
mentawai = ggplot() + 
  geom_raster(data = interpolated_df_mentawai, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_mentawai_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_mentawai, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
mentawai

#=== Pesisir Selatan
kec_pesisir_selatan<-kec_cropped[kec_cropped$KECAMATAN %in%c("IV JURAI","IV NAGARI BAYANG UTARA","AIRPURA","BASA AMPEK BALAI TAPAN","BATANG KAPAS","BAYANG","KOTO XI TARUSAN","LINGGO SARI BAGANTI","LENGAYANG","LUNANG","PANCUNG SOAL","RANAH AMPEK HULU TAPAN","RANAH PESISIR","SILAUT","SUTERA"),]
kec_pesisir_selatan_sf <- st_as_sf(kec_pesisir_selatan)
kec_pesisir_selatan_sf_utm1 <- st_transform(kec_pesisir_selatan_sf, crs = 32648)

Earthquakes_Pred_Pesisir_Selatan <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_pesisir_selatan_sf_utm1)$xmin) & x <= max(st_bbox(kec_pesisir_selatan_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_pesisir_selatan_sf_utm1)$ymin) & y <= max(st_bbox(kec_pesisir_selatan_sf_utm1)$ymax))

Earthquakes_Pred_Pesisir_Selatan$x <- jitter(Earthquakes_Pred_Pesisir_Selatan$x, amount = 0.0001)
Earthquakes_Pred_Pesisir_Selatan$y <- jitter(Earthquakes_Pred_Pesisir_Selatan$y, amount = 0.0001)

interpolated_PS <- with(Earthquakes_Pred_Pesisir_Selatan, 
                        interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_PS <- as.data.frame(expand.grid(x = interpolated_PS$x, y = interpolated_PS$y))
interpolated_df_PS$pred_mean <- as.vector(interpolated_PS$z)
interpolated_df_PS <- na.omit(interpolated_df_PS)  # Menghapus nilai NA

interpolated_df_PS$title = "Pesisir Selatan"
ps = ggplot() + 
  geom_raster(data = interpolated_df_PS, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_pesisir_selatan_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_pesisir_selatan, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
ps

#=== Padang lawas and Mandailing natal
kab_sumut_wgs84 <- st_transform(kab_sumut, crs = st_crs(kec))
bbox_kab_sumut <- st_bbox(kab_sumut_wgs84)
kec_cropped2 <- st_crop(kec_valid, bbox_kab_sumut)

kec_plm<-kec_cropped2[kec_cropped2$KECAMATAN %in%c("AEK NABARA BARUMUN","BARUMUN","BARUMUN BARAT","BARUMUN BARU","BARUMUN SELATAN","BARUMUN TENGAH","BATANG LUBU SUTAM","HURISTAK","HUTA RAJA TINGGI","LUBUK BARUMUN","SIHAPAS BARUMUN","SOSA","SOSA JULU","SOSA TIMUR","SOSOPAN","ULU BARUMUN","ULU SOSA",
                                                   "BATAHAN","BATANG NATAL","BUKIT MALINTANG","HUTA BARGOT","KOTANOPAN","LEMBAH SORIK MARAPI","LINGGA BAYU","MUARA BATANG GADIS","MUARA SIPONGI","NAGA JUANG","NATAL","PEKANTAN","PENYABUNGAN BARAT","PENYABUNGAN KOTA","PENYABUNGAN SELATAN","PENYABUNGAN TIMUR","PENYABUNGAN UTARA","PUNCAK SORIK MARAPI","RANTO BAEK","SIABU","SINUNUKAN","TAMBANGAN","ULU PUNGKUT"),]
kec_plm_sf <- st_as_sf(kec_plm)
kec_plm_sf_utm1 <- st_transform(kec_plm_sf, crs = 32648)

Earthquakes_Pred_PLM <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_plm_sf_utm1)$xmin) & x <= max(st_bbox(kec_plm_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_plm_sf_utm1)$ymin) & y <= max(st_bbox(kec_plm_sf_utm1)$ymax))

Earthquakes_Pred_PLM$x <- jitter(Earthquakes_Pred_PLM$x, amount = 0.0001)
Earthquakes_Pred_PLM$y <- jitter(Earthquakes_Pred_PLM$y, amount = 0.0001)

interpolated_PLM <- with(Earthquakes_Pred_PLM, 
                         interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_PLM <- as.data.frame(expand.grid(x = interpolated_PLM$x, y = interpolated_PLM$y))
interpolated_df_PLM$pred_mean <- as.vector(interpolated_PLM$z)
interpolated_df_PLM <- na.omit(interpolated_df_PLM)  # Menghapus nilai NA

interpolated_df_PLM$title = "Padang Lawas & Mandailing Natal"
plm = ggplot() + 
  geom_raster(data = interpolated_df_PLM, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_plm_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_plm, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
plm

#=== Nias 
kec_nias<-kec_cropped2[kec_cropped2$KECAMATAN %in%c("AFULU","ALASA","ALASA TALUMUZOI","LAHEWA","LAHEWA TIMUR","LOTU","NAMOHALU","SAWO","SITOLU ORI","TUGALA OYO","TUHEMBERUA","GUNUNGSITOLI","GUNUNGSITOLI ALO'OA","GUNUNGSITOLI BARAT","GUNUNGSITOLI IDANOI","GUNUNGSITOLI SELATAN","GUNUNGSITOLI UTARA","BAWOLATO","BOTOMUZOI","GIDO","HILIDUHO","HILISERANGKAI","IDANOGAWO","MA'U","SOGAE'ADU","SOMOLO-MOLO","ULUGAWO","LAHOMI","LOLOFITU MOI", "MANDREHE","MANDREHE BARAT","MANDREHE UTARA","MORO'O","SORIMBU","UKU MORO'O"),]
kec_nias_sf <- st_as_sf(kec_nias)
kec_nias_sf_utm1 <- st_transform(kec_nias_sf, crs = 32648)

Earthquakes_Pred_Nias <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_nias_sf_utm1)$xmin) & x <= max(st_bbox(kec_nias_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_nias_sf_utm1)$ymin) & y <= max(st_bbox(kec_nias_sf_utm1)$ymax))

Earthquakes_Pred_Nias$x <- jitter(Earthquakes_Pred_Nias$x, amount = 0.0001)
Earthquakes_Pred_Nias$y <- jitter(Earthquakes_Pred_Nias$y, amount = 0.0001)

interpolated_Nias <- with(Earthquakes_Pred_Nias, 
                          interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_Nias <- as.data.frame(expand.grid(x = interpolated_Nias$x, y = interpolated_Nias$y))
interpolated_df_Nias$pred_mean <- as.vector(interpolated_Nias$z)
interpolated_df_Nias<- na.omit(interpolated_df_Nias)  # Menghapus nilai NA

interpolated_df_Nias$title = "Nias"
nias = ggplot() + 
  geom_raster(data = interpolated_df_Nias, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_nias_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_nias, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
nias

#=== South Nias
kec_nias_sel<-kec_cropped2[kec_cropped2$KECAMATAN %in%c("AMANDRAYA","ARAMO","BORONADU","FANAYAMA","GOMO","HIBALA","HILIMEGAI","HILISALAWA'AHE","HURUNA","IDANOTAE","LAHUSA","LOLOMATUA","LOLOWAU","LUAHAGUNDRE MANIAMOLO","MANIAMOLO","MAZINO","MAZO","O'O'U","ONOHAZUMBA","ONOLALU","PULAU-PULAU BATU","PULAU-PULAU BARAT","PULAU-PULAU BATU TIMUR", "PULAU-PULAU BATU UTARA","SIDUA'ORI","SIMUK",
                                                        "SOMAMBAWA","SUSUA","TANAH MASA","TOMA","ULUNOYO", "ULU IDANOTAE","UMBUNASI","ULUSUSUA"),]
kec_niassel_sf <- st_as_sf(kec_nias_sel)
kec_niassel_sf_utm1 <- st_transform(kec_niassel_sf, crs = 32648)

Earthquakes_Pred_Nias_Sel <- Earthquakes_Pred_dec %>%
  filter(x >= min(st_bbox(kec_niassel_sf_utm1)$xmin) & x <= max(st_bbox(kec_niassel_sf_utm1)$xmax) &
           y >= min(st_bbox(kec_niassel_sf_utm1)$ymin) & y <= max(st_bbox(kec_niassel_sf_utm1)$ymax))

Earthquakes_Pred_Nias_Sel$x <- jitter(Earthquakes_Pred_Nias_Sel$x, amount = 0.0001)
Earthquakes_Pred_Nias_Sel$y <- jitter(Earthquakes_Pred_Nias_Sel$y, amount = 0.0001)

interpolated_Nias_Sel <- with(Earthquakes_Pred_Nias_Sel, 
                              interp(x = x, y = y, z = pred_mean, nx = 200, ny = 200))

interpolated_df_Nias_Sel <- as.data.frame(expand.grid(x = interpolated_Nias_Sel$x, y = interpolated_Nias_Sel$y))
interpolated_df_Nias_Sel$pred_mean <- as.vector(interpolated_Nias_Sel$z)
interpolated_df_Nias_Sel<- na.omit(interpolated_df_Nias_Sel)  # Menghapus nilai NA

interpolated_df_Nias_Sel$title = "Nias Selatan"
nias_sel = ggplot() + 
  geom_raster(data = interpolated_df_Nias_Sel, aes(x = x, y = y, fill = pred_mean), interpolate = TRUE) + 
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"), limits = c(0, max_pred)) + 
  geom_sf(data = kec_niassel_sf_utm1, fill = NA, color = "black", linewidth = 0.3) + 
  geom_sf_text(data = kec_nias_sel, aes(label = KECAMATAN), size = 1.5, color = "black", fontface = "bold") + 
  coord_sf(crs = 32648) + 
  theme_bw() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.7),
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
  ) + 
  labs(
    x = "Longitude", y = "Latitude", fill = "Earthquake Magnitude"
  ) +
  facet_grid(. ~ title)
nias_sel

#================== 9.3 High-Resolution of Exceedance Probability
#================9.3.1 Overall: By Month
ggplot() + 
  geom_tile(data = Earthquakes_Pred, aes(x = x, y = y, fill = pred_prob)) + 
  geom_contour(data = Earthquakes_Pred, aes(x = x, y = y, z = pred_mean, colour = after_stat(level)), 
               linewidth = 0.2, colour = "red", alpha = 0.5) +
  scale_fill_gradientn(colours = c("gray99", "yellow", "red")) + 
  facet_wrap(~ Time, labeller = as_labeller(Month_names), ncol = 4) + 
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  geom_sf(data = Sumatra_sf_utm1, fill = NA, size = 0.01) +  
  theme(legend.position = "bottom", text = element_text(size = 19)) + 
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1, 
                               frame.colour = "black",
                               barwidth = 20,
                               barheight = 1.5)) +
  theme(axis.text.x = element_blank(), #remove x axis labels
        axis.ticks.x = element_blank(), #remove x axis ticks
        axis.text.y = element_blank(),  #remove y axis labels
        axis.ticks.y = element_blank()  #remove y axis ticks
  ) + 
  labs(fill = "Earthquakes \nExceedance Probability") + 
  ggtitle("Earthquakes Exceedance Probability")

#================9.3.2 Overall: By Province
max_prob = max(Earthquakes_Pred$pred_prob); max_prob
ggplot() +    
  geom_tile(data = Earthquakes_Pred_dec, aes(x = x, y = y, fill = pred_prob)) +    
  geom_contour(data = Earthquakes_Pred_dec, aes(x = x, y = y, z = pred_mean, colour = stat(level)), 
               size = 0.2, colour = "red", alpha = 0.5) +   
  scale_fill_gradientn(colours = c("gray99", "yellow", "red"),limits = c(0, max_prob)) +   
  theme_bw() +    
  ylab("Latitude") +    
  xlab("Longitude") +    
  geom_sf(data = Sumatra_sf_utm1, fill = NA, size = 0.01) +     
  geom_sf_text(data =  Sumatra_sf_utm1, aes(label = NAME_1), size = 3, color = "black", fontface = "bold") +
  theme(legend.position = "bottom", 
        text = element_text(size = 12)) +   
  theme(
    panel.background = element_rect(fill = "white", color = NA), 
    panel.grid.major = element_line(color = "grey95", linetype = "solid"),
    panel.grid.minor = element_line(color = "grey95", linetype = "solid"),
    plot.title = element_text(hjust = 0.5, size = 12),  
    panel.grid = element_blank()                        
  )+
  guides(fill = guide_colorbar(title.position = "left", title.vjust = 1, 
                               frame.colour = "black", 
                               barwidth = 20, 
                               barheight = 1.5))+
  facet_grid(. ~ title)