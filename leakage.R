#### Setup ####
library(bife)
library(plyr)
library(tidyverse)
library(survival)
library(jsonlite)
library(texreg)
library(MASS)
library(stats)
library(matrixStats)
library(magrittr)
library(reshape2)
library(sf)
library(lubridate)
library(rnaturalearth)
library(rnaturalearthdata)
library(scales)
library(rgdal)
library(kableExtra)
library(zoo)
library(ggsn)
library(cowplot)
library(raster)
library(ggnewscale)
library(janitor)
library(ggridges)
library(rgeos)
install.packages("rnaturalearthhires", repos = "http://packages.ropensci.org", type = "source")

select <- dplyr::select

#### Load/clean data ####
## Load primary datasets
full.df <- read.csv("input/long_kali.csv")
mill_df <- read_csv("input/master_mill_data.csv")

## Define parameter for conversion of points (2x2km sample grid) to area (4km2 representation)
area_wght <- 4

## Filter to points starting as forest (drops from 106889 to 59513 points)
df <- full.df %>%
  filter(for_2000 != 0)

## Clean up variables
df$defor <- !df$forest
df$pid <- as.numeric(as.factor(df$propid))

## Define kh to be used
df$kh <- df$kh_2005

## ID points belonging to RSPO member company with some certification elsewhere
rspo_member <- df %>%
  group_by(grp_id) %>%
  summarize(on_rspo = max(g_cert_shr) > 0)
df <- full_join(df, rspo_member, by = "grp_id")

## Create op_type variable to indicate whether points are not on oil palm concession (1),
## on an rspo member concession (2), or on an ever-certified concession
df$on_op <- df$conc_class==3
df$op_type = 1
df <- df %>%
  mutate(op_type = replace(op_type, on_op==1, 2),
         op_type = replace(op_type, cert==1, 3))

## Filter sample to management units with some actual variation in deforestation
## These are the only units that bife approach can use. Drops sample further to 48157 samples
n_outcomes = df  %>%
  filter(!is.na(defor)) %>% 
  group_by(pid) %>%     
  summarise(n_out = n_distinct(defor))
keep <- as.list(n_outcomes[(n_outcomes$n_out >1), "pid"])$pid
sample <- df %>% filter(pid %in% keep)

## Reorder sample to match what comes out of BIFE
sample <- sample[order(sample$pid, sample$year),]

#### Regressions ####
## Define helper function to run bife
bc_APE <- function(bife_model){
  apes <- get_APEs(bias_corr(bife_model))
  apes$nobs <- nobs
  apes$npts <- npts
  return(apes)
}

## Descriptive results - likelihood of deforestation on different land use zones (T2)
npts <- sample$sid %>% unique() %>% length()
nobs <- length(sample$defor[!is.na(sample$defor)])

print(npts)
print(nobs)

## Main regression results - building up to full model (T3)
cert_ape <- bife(defor ~ cert_now + 
                   factor(kh) + factor(year) | pid, data = sample)  %>% bc_APE()

certkh_ape <- bife(defor ~ cert_now:factor(kh) + 
                   factor(kh) + factor(year) | pid, data = sample)  %>% bc_APE()

pcc_ape <- bife(defor ~ cert_now + pccp + 
                      factor(kh) + factor(year) | pid, data = sample) %>% bc_APE()

pcckh_ape <- bife(defor ~ cert_now:factor(kh) + pccp:factor(kh) + 
                        factor(kh) + factor(year) | pid, data = sample) %>% bc_APE()

gc_ape <- bife(defor ~ cert_now + g_cert_shr + 
                      factor(kh) + factor(year) | pid, data = sample) %>% bc_APE()

gckh_ape <- bife(defor ~ cert_now:factor(kh) + g_cert_shr:factor(kh) + 
                        factor(kh) + factor(year) | pid, data = sample) %>% bc_APE()

all_mod <- bife(defor ~ cert_now + pccp + g_cert_shr + 
                      factor(kh) + factor(year) | pid, data = sample)
all_bc <- bias_corr(all_mod)
all_ape <- all_mod %>% bc_APE()

allkh_mod <- bife(defor ~ cert_now:factor(kh) + pccp:factor(kh) + g_cert_shr:factor(kh) + 
                    factor(kh) + factor(year) | pid, data = sample)
allkh_bc <- bias_corr(allkh_mod)
allkh_ape <- allkh_mod %>% bc_APE()

#### Simulation ####
## Define helper functions for simulation
point_persist <- function(simulation) {
  persist_prob <- 1 - simulation
  cumprod <- persist_prob %>%
    group_by(sample$sid) %>%
    summarise_all(list(cum_defor = prod)) %>%
    rename(sid = "sample$sid")
  return(cumprod)
}

agg_sims <- function(simulated, sample, sum_cats) {
  cumprod <- point_persist(simulated)
  attributes <- sample %>%
    select(sid, sum_cats) %>%
    unique()
  
  defor_df <- full_join(cumprod, attributes, by = "sid")
  defor_df <- defor_df %>% drop_na()
  
  sim_agg <- defor_df %>%
    group_by_at(sum_cats) %>%
    summarise_at(vars(matches("cum_defor")), list(defor = sum)) %>%
    ungroup()
  
  return(sim_agg)
}

smry_stats <- function(agg.sim) {
  simtab <- agg.sim %>%
    select(contains("_cum_defor")) * area_wght
  agg.sim$mean <- simtab %>%
    rowMeans()
  agg.sim$std <- simtab %>%
    as.matrix %>%
    rowSds()
  stats_smry <- agg.sim %>%
    select(sum_cats, mean, std)
  return(stats_smry)  
}

## Define regression to be used for simulations
sim_mod <- allkh_bc

## Attach alphas to sample
alpha <- sim_mod$alpha 
pids <- names(alpha)
alpha <- alpha %>% as.numeric()
alpha <- tibble(pid = pids, alpha = alpha)
alpha$pid <- alpha$pid %>% as.numeric()
sample <- full_join(sample, alpha, by = "pid")

## Prepare design matrix for baseline simulation
formula <- sim_mod$formula
X.bl <- model.matrix(formula, sample)[,-1]

## Prepare design matrix for counterfactual simulation (no certification)
data.cf <- sample
data.cf$cert_now = 0
data.cf$pccp = 0
data.cf$g_cert_shr = 0
X.cf <- model.matrix(formula, data.cf)[,-1]

## Draw parameters for simulation
mod_adj <- sim_mod
sum_cats <- c("op_type", "kh")
beta <- sim_mod$coefficients
beta_cov <- vcov(sim_mod)
n <- 1000
ptm <- proc.time()
beta_draw <- t(mvrnorm(n, mu = beta, Sigma = beta_cov))

## Run simulations
simulated.bl <- data.frame(matrix(ncol = n, nrow = dim(sample)[1]))
simulated.cf <- data.frame(matrix(ncol = n, nrow = dim(sample)[1]))
for (i in 1:n) {
  mod_adj$coefficients <- matrix(beta_draw[,i])
  simulated.bl[,i] <- predict(mod_adj, X_new = X.bl, alpha_new = sample$alpha, type = "response")
  simulated.cf[,i] <- predict(mod_adj, X_new = X.cf, alpha_new = sample$alpha, type = "response")
  }

## Calculate deforestation rates
agg.cf <- agg_sims(simulated.cf, sample, sum_cats)
smry.cf <- smry_stats(agg.cf)
agg.bl <- agg_sims(simulated.bl, sample, sum_cats)
smry.bl <- smry_stats(agg.bl)
agg.dif <- agg.bl - agg.cf
agg.dif[,1:length(sum_cats)] = agg.cf[,1:length(sum_cats)]
smry.dif <- smry_stats(agg.dif)
proc.time() - ptm

start_for <- sample %>%
  filter(year==2004) %>%
  group_by_at(sum_cats) %>%
  summarise(for_area = sum(forest, na.rm = TRUE) * area_wght)
agg.shr <- agg.dif
agg.shr[, grep("_cum_defor", names(agg.shr))] <- agg.shr[, grep("_cum_defor", names(agg.shr))] / start_for$for_area * 100


#### Tables ####
## Extract value from 'bife' models
extract.tex <- function(bife.ape) {
  sum.ape <- summary(bife.ape)
  coef.names <- row.names(sum.ape)
  coef <- sum.ape[,1]
  se <- sum.ape[,2]
  pvalues <- sum.ape[,4]
  model.name <- ""
  nobs <- nobs
  npts <- npts
  gof <- c(nobs, npts)
  gof.names <- c("N. observations", "N. points") 
  gof.decimal <- c(FALSE,FALSE)
  tex <- createTexreg(coef.names, coef, se, pvalues, gof.names = gof.names, gof = gof,gof.decimal = gof.decimal)
  return(tex)
}

## Table 2: Main spillover regressions
cert.bife <- extract.tex(cert_ape)
certkh.bife <- extract.tex(certkh_ape)
pcc.bife <- extract.tex(pcc_ape)
pcckh.bife <- extract.tex(pcckh_ape) 
gc.bife <- extract.tex(gc_ape) 
gckh.bife <- extract.tex(gckh_ape) 
all.bife <- extract.tex(all_ape) 
allkh.bife <- extract.tex(allkh_ape) 

# Create regression table
map <- list("cert_now" = "Certified",
            "cert_now:factor(kh)0" = "Certified x Not forest estate",
            "cert_now:factor(kh)1" = "Certified x Forest estate",
            "g_cert_shr" = "Certified share of parent company's holdings",
            "g_cert_shr:factor(kh)FALSE" = "Certified share of parent company's holdings x Not forest estate",
            "factor(kh)0:g_cert_shr" = "Certified share of parent company's holdings x Not forest estate",
            "g_cert_shr:factor(kh)TRUE" = "Certified share of parent company's holdings x Forest estate",
            "factor(kh)1:g_cert_shr" = "Certified share of parent company's holdings x Forest estate",
            "pccp" = "Certified share of local mill capacity",
            "factor(kh)0:pccp" = "Certified share of local mill capacity x Not forest estate",
            "factor(kh)1:pccp" = "Certified share of local mill capacity x Forest estate")

htmlreg(list(cert.bife,certkh.bife, gc.bife, gckh.bife, pcc.bife, pcckh.bife, all.bife, allkh.bife), 
       custom.coef.map = map,use.packages = FALSE,scalebox = 0.65,label = "tab3", stars = c(0.01, 0.05, 0.1),
       digits =4, caption="Spillovers",caption.above=TRUE,no.margin = TRUE,column.spacing=0,
       file = "output/tbl2_spillover_regression.doc")

## Table S2 - Areas of forests and deforestation
# Calculating forest areas within and outside certified supply shed area
cert_st <- sample %>%
  select(sid,pccp) %>%
  group_by(sid) %>%
  summarize(in_ss = max(pccp) > 0)

lu_areas <-  sample %>%
  filter(year==2004, !is.na(defor)) %>%
  right_join(cert_st,by="sid") %>%
  filter(op_type==1 | op_type==2 | op_type ==3) %>%
  select(sid,op_type,kh_2005,kh_2018,in_ss) %>%
  group_by(in_ss,op_type,kh_2005,kh_2018) %>%
  summarise(area = n()*area_wght/1000) %>%
  ungroup()

# Reformat to create table
lu_areas_tbl <- lu_areas %>%
  add_column(conc_type="") %>%
  mutate(conc_type = ifelse(op_type ==1,"Outside oil palm concession",
                              ifelse(op_type==2,"Non certified oil palm concession",
                                     ifelse(op_type==3,"Certified oil palm concession",conc_type)))) %>%
  mutate(in_ss = ifelse(in_ss == "TRUE","Inside 2016 certified supply shed","Outside certified supply shed")) %>%
  mutate(kh_2005 = ifelse(kh_2005==0,"APL","Forest Estate")) %>%
  mutate(kh_2018 = ifelse(kh_2018==0,"APL","Forest Estate")) %>%
  mutate(zone_change =   paste(kh_2005,kh_2018, sep = ' to ')) %>%
  select(in_ss,conc_type,zone_change,area) %>%
  spread(zone_change,area) %>%
  select(in_ss,conc_type,`Forest Estate to Forest Estate`,`Forest Estate to APL`,`APL to APL`,`APL to Forest Estate`) %>%
  mutate(Total = `Forest Estate to Forest Estate`+ `Forest Estate to APL` +`APL to APL`+`APL to Forest Estate`) %>%
  drop_na() 

# Creating table using kable (for forest areas inside and outside 2016 certified supply shed)
tbl3<- kable(lu_areas_tbl[c(1,2,3,4,5),2:7], "html",align = "r",caption=NA,row.names=FALSE,col.names=c("","","","","",""),booktabs=TRUE) %>%
  #kable_styling(position = "center",full_width = FALSE,latex_options = "hold_position")%>%
  pack_rows(index = c("Inside 2016 certified supply shed" = 3, "Outside 2016 certified supply shed" = 2),indent = F) %>%
  add_header_above(c("", paste0("(thousand km", "<sup>", 2, "</sup>",")"),paste0("(thousand km", "<sup>", 2, "</sup>",")"),
                     paste0("(thousand km", "<sup>", 2, "</sup>",")"),paste0("(thousand km", "<sup>", 2, "</sup>",")"),
                     paste0("(thousand km", "<sup>", 2, "</sup>",")")),bold = T,align="c",escape = FALSE) %>%
  add_header_above(c("Concession type","Forest Estate to Forest Estate","Forest Estate to APL","APL to APL","APL to Forest Estate","Total"),bold=T,align="c",escape = FALSE, line = F)

write_file(tbl3,"output/tbl3_zone_change_areas_ss.doc")

#### Figures ####

## Figure 3 ##
# Reformatting data for figure
certyr_df <- mill_df %>%
   select(db_id,ci_year) %>%
   mutate(ci_year = ci_year + 1)

# Get production (capacity) data in long format
prodlong_df <- mill_df %>%
  select(db_id,`2004`:`2016`) %>% 
  gather(year,cap,-db_id)

certyr_long_df <- mill_df %>%
  select(db_id,ci_year) %>%
  left_join(prodlong_df, by = "db_id", all.x = TRUE) %>%
  group_by(db_id) %>%
  mutate(cert = na.locf(ci_year, na.rm = F),ci_year = ci_year + 1) 

certyr_long_df$cert <- ifelse(certyr_long_df$year >= certyr_long_df$ci_year,1,0)

prodmerge_df = merge(prodlong_df,certyr_long_df,by=c("db_id","year"),all.y=TRUE)
prodmerge_df[is.na(prodmerge_df)] <- 0

certst_df = merge(mill_df,prodlong_df,by="db_id",all.y=TRUE)
certst_df$cert_status[is.na(certst_df$cert_status)] <- 0
certst_df$prov=separate(certst_df ,db_id, sep = "\\-\\-*",c("country_id","prov"))[,2]

# Only select data for Kalimantan
certst_df <- certst_df[grep("K", certst_df$prov),]

# Convert dates from year to ymd (year-month-day)
certst_df$year <- ymd(certst_df$year, truncated = 2L)

certmem_df <- mill_df %>% 
  select(grp_id,rspo_member,year = year_join_rspo) %>%
  drop_na() %>%
  distinct(grp_id, .keep_all= TRUE) %>%
  mutate(year = ymd(year, truncated = 2L))

certst_df$cert_status_2 = ifelse(certst_df$cert_status == 'not-certified' & certst_df$grp_id %in% unique(certmem_df$grp_id) ,'non-certified member',certst_df$cert_status)
nc_mem_df <- subset(certst_df,cert_status_2 =='non-certified member')

# Calculate total capacity of non-certified members
nc_mem_df = merge(certst_df,certmem_df,by="grp_id",all.X=TRUE)
nc_mem_df$cert_status_2 = ifelse(nc_mem_df$year.y > nc_mem_df$year.x, 'certified member',nc_mem_df$cert_status_2)

ncmem_df <- nc_mem_df %>%
  group_by(prov,year.x,cert_status_2) %>%
  summarize(total_caps = sum(cap,na.rm=T))

ncmem_df <- subset(ncmem_df,cert_status_2 =='non-certified member')

noncertmem_df <- ncmem_df %>%
  select(prov=prov,year=year.x,cert_status=cert_status_2,total_caps)

certcap_df = certst_df %>%
        group_by(prov,year,cert_status) %>%
        summarize(total_caps = sum(cap,na.rm=T))

certpot_df = subset(certcap_df,cert_status == "certified")

# Extract province code
prodmerge_df$prov=separate(prodmerge_df ,db_id, sep = "\\-\\-*",c("country_id","prov_code"))[,2]
prodmerge_df <- prodmerge_df[grep("K", prodmerge_df$prov),]

# Get total mill capacity by province and year
mill_cap_df = prodmerge_df %>%
          group_by(prov,year,cert) %>%
          summarize (total_caps = sum(cap.x))

mill_cap_df$cert_status = ifelse(mill_cap_df$cert == 1, 'Certified','Not-certified')
mill_cap_df$year <- ymd(mill_cap_df$year, truncated = 2L)

certnc_df <- mill_cap_df %>%
  select(prov,year,cert_status,total_caps)

# Merge data frames
prov_caps_df = rbind(noncertmem_df,certnc_df)

# Convert from long to wide
prov_caps_wide_df <- dcast(prov_caps_df, prov + year ~ cert_status, fill=0)
prov_caps_wide_df$`Not-certified` = prov_caps_wide_df$`Not-certified` - prov_caps_wide_df$`non-certified member`

# Rename provinces
prov_caps_wide_df$prov[prov_caps_wide_df$prov == "KB"] <- "West Kalimantan"
prov_caps_wide_df$prov[prov_caps_wide_df$prov == "KT"] <- "Central Kalimantan"
prov_caps_wide_df$prov[prov_caps_wide_df$prov == "KS"] <- "South Kalimantan"
prov_caps_wide_df$prov[prov_caps_wide_df$prov == "KI"] <- "East Kalimantan"
prov_caps_wide_df$prov[prov_caps_wide_df$prov == "KU"] <- "North Kalimantan"

# Convert from wide to long
prov_caps_long_df <- melt(prov_caps_wide_df, id.vars = c("prov", "year"), variable.name = "cert")
prov_caps_long_df$cert <- as.character(prov_caps_long_df$cert)
prov_caps_long_df$year <- as.Date(prov_caps_long_df$year,format="%Y-%M-%D")

# Renaming cert status to get sorting for plot
prov_caps_long_df$cert[prov_caps_long_df$cert == "Certified"] <- "1 - RSPO certified"
prov_caps_long_df$cert[prov_caps_long_df$cert == "non-certified member"] <- "2 - RSPO member but not certified"
prov_caps_long_df$cert[prov_caps_long_df$cert == "Not-certified"] <- "3 - Not covered by RSPO"

# Stacked Area Plot
cs_plot <- ggplot(prov_caps_long_df, aes(x=year, y=value, fill=cert))+
           geom_area(stat='identity', position='stack',colour=NA)+ xlab("\n") +
           ylab("Total mill capacity (Tonnes FFB / hour)") +
           facet_wrap("prov",nrow=1,scales='free_x')+
           scale_x_date(expand=c(0,0), breaks=date_breaks('2 years'), labels=date_format('%Y'))+
           scale_y_continuous(expand=c(0,1),limits=c(0,6500),breaks=c(0,1000,2000,3000,4000,5000,6000,7000))+
           scale_fill_manual(values=c('darkturquoise','deeppink','gold'),name="Certification Status",labels=c('RSPO Certified','RSPO Member but not certified','Not covered by RSPO')) +
           theme_bw()+
           guides(fill = guide_legend(title.position = "top",title.hjust =0.5))+
           theme(axis.text.x=element_text(vjust=0,hjust=0.1,size=8,angle = 45),
           axis.title.y=element_text(size=8.5),
           panel.spacing = unit(0.3, "in"),
           panel.grid.major = element_blank(),
           panel.grid.minor = element_blank(),
           plot.margin=unit(c(0.2,0.2,0.2,0.2),"in"),
           strip.background=element_rect(fill=NA, color=NA),
           legend.position="bottom",
           legend.direction='horizontal',
           legend.text = element_text(size = 8),
           legend.title = element_text(size = 9),
           legend.justification="center",
           legend.margin=margin(0,0,0,0),
           legend.box.margin=margin(3,3,3,3))
cs_plot

# save plot to png
ggsave(cs_plot,file="output/fig03_cert_status.png",dpi=500,w=9,h=5,unit="in",type="cairo-png")

#### Figure 4 
## Impact of certification - simulated densities
type_totals <- agg.dif %>%
  group_by(op_type) %>%
  summarize_all(sum)
type_totals$kh <- -1
net_change <- agg.dif %>% 
  summarize_all(sum)
net_change$op_type <- 0
net_change$kh <- -1
kh_totals <- agg.dif %>%
  group_by(kh) %>%
  summarize_all(sum)
kh_totals$op_type <- 0

all.dif <- rbind(agg.dif, net_change, type_totals, kh_totals)

all.dif %>% 
  filter(op_type==0, kh==1) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE) %>% 
  quantile(c(0.025, 0.5, 0.975))

long.dif <- all.dif %>% 
  melt(id.vars = sum_cats, value.name = "dif") %>%
  unite(use_class, sum_cats, sep = '_', remove = FALSE)  %>%
  mutate(kh = as.numeric(kh))

long.dif <- long.dif %>%
  unite(use_class, c("op_type"), sep = '_', remove = FALSE)

## Box plot
ioc_plot <- ggplot(long.dif, aes(x = use_class, y = dif, fill = as.character(kh))) + 
            geom_boxplot(aes(middle = mean(dif)),notch=FALSE) +
            coord_flip() +
            geom_vline(xintercept = 1.5, linetype = "solid", color = "black") +
            scale_point_color_hue(l = 40) +
            scale_discrete_manual(aesthetics = "point_fill", values = c(0, 1)) +
            ylab("Change in forest area"~(km^2)) + xlab("Land use category\n")+
            geom_hline(yintercept = 0, color = "black", size = 1,linetype="longdash") +
            scale_fill_cyclical(values = c("grey40", "deeppink", "darkturquoise"),
            name = "Land use zone", guide = "legend",
            labels = c( "Sum", "APL", "Forest estate")) +
            theme_bw()+
            scale_x_discrete(breaks=c(0, 1, 2, 3),
            labels=c("Total deforestation\nacross landscape",
                     "Outside oil palm\nconcessions", 
                     "Never certified\noil palm concessions",
                     "Certified oil palm\nconcessions")) +
            theme(legend.position="bottom",
            legend.box="horizontal",
            legend.background = element_blank(),
            legend.box.background = element_blank(),
            axis.text.x=element_text(size=9,face="bold"),
            axis.text.y=element_text(size=9,face="bold"),
            plot.title = element_text(hjust = 0.5)) +
            guides(fill=guide_legend(title.position="top", 
            title.hjust =0.5))

ioc_plot

# save plot to png
ggsave(ioc_plot,file="output/fig04_impactofdefor_boxplot.png",dpi=500,w=9,h=5,unit="in",type="cairo-png")


## Figure 5
persist.cf <- point_persist(simulated.cf)
persist.bl <- point_persist(simulated.bl)
sim_cols <- grep("cum_defor", names(persist.bl))
defor.dif <- persist.bl
defor.dif[,sim_cols] <- (1-persist.bl[,sim_cols]) - (1-persist.cf[,sim_cols])
defor.dif$defor_dif <- defor.dif[,sim_cols] %>%
  rowMeans() * 100
defor.dif <- defor.dif %>%
  select(sid, defor_dif)

sample_sf <-  sample %>% 
  filter(year == 2004) %>%
  select(sid, lat, lon) %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326)

sample_sf <- full_join(sample_sf, defor.dif, by = "sid") %>%
  drop_na()

coast_sf <- ne_coastline(scale = "large", returnclass = "sf")
countries_sf <- ne_countries(scale = "large", returnclass = "sf")
ocean <- ne_download(scale = 10, type = 'ocean', category = 'physical')
ocean_sf <- st_as_sf(ocean)

# reproject to equal area projection
indonesian_crs <- "+proj=cea +lon_0=115.0 +lat_ts=0 +x_0=0 +y_0=0 +ellps=WGS84 +units=m +no_defs"

sample_proj = st_transform(sample_sf,indonesian_crs)
raster_template = raster(extent(sample_proj), resolution = 2000,
                         crs = st_crs(sample_proj)$proj4string)

# rasterize samples
sample_ras = rasterize(sample_proj, raster_template, field = "defor_dif")
sr <- "+proj=longlat +datum=WGS84 +no_defs" 

# reproject back to wgs84 for map
sample_ras <- projectRaster(sample_ras, crs = sr)
sample_spdf <- as(sample_ras, "SpatialPixelsDataFrame")
sample_df <- as.data.frame(sample_spdf)
colnames(sample_df) <- c("value", "x", "y")

# create class cuts for map and labels
options(scipen = 999) # ensure numbers are not in scientific format
cuts <- c(-35, -15, -7, -3.5, -1.25, -0.0001, 0.0001, 1.25, 3.5, 7, 15, 35)
no_classes <- 11
labels <- c()

labels <- c()
for(idx in 1:length(cuts)){
  labels <- c(labels, paste0(round(cuts[idx],8), 
                             " to ",  round(cuts[idx + 1], 5)))
}

# remove the last label 
# because that would be something like "35 - NA"
labels <- labels[1:length(labels)-1]

# create dataframe with deforestation probability classes
sample_df$value <- cut(sample_df$value,
                  breaks = cuts, 
                  labels = labels, 
                  include.lowest = T)

# create more legible legend of probabilities
sample_df <- sample_df %>%
  mutate(label_value = as.character(as.factor(value))) %>%
  mutate(label_value = ifelse(label_value == "-1.25 to -0.0001","-1.25 to 0",label_value),
  label_value = ifelse(label_value == "-0.0001 to 0.0001","0",label_value),
  label_value = ifelse(label_value == "0.0001 to 1.25","0 - 1.25",label_value)) %>%
  mutate(label_value = as.factor(label_value))

# detach raster package to use ggsn for scale bar on maps
detach("package:raster", unload=TRUE)

# Set up theme
theme_map <- theme(text = element_text(family = "Arial", color = "black"),
                   panel.background = element_rect(fill="transparent"),
                   panel.border = element_rect(colour = "black",fill=F,size=1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   plot.title = element_text(hjust = 0.5),
                   axis.line = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank(),
                   axis.title.x = element_blank(),
                   axis.title.y = element_blank(),
                   legend.position = c(0.2, 0.78),
                   legend.direction = 'vertical',
                   legend.key.height = unit(6, "pt"),
                   legend.key.width = unit(10, "pt"),
                   legend.text = element_text(size = 8.5),
                   legend.title = element_text(size = 9.5),
                   legend.spacing.x = unit(0.1,'cm'),
                   legend.title.align = 0,
                   legend.text.align = 0,
                   legend.background = element_rect(colour=NA,fill=NA))

# reordering map labels
sample_df$label_value <- ordered(sample_df$label_value, levels = c("-35 to -15","-15 to -7","-7 to -3.5",
                          "-3.5 to -1.25","-1.25 to 0","0","0 to 1.25",
                          "1.25 to 3.5","3.5 to 7","7 to 15","15 to 35"))

# create overview map A
oa_map <- ggplot(sample_df) +  
  geom_tile(aes(x=x, y=y, fill=label_value,color=label_value),alpha=0.8,size=0.0001) + 
  coord_sf() +
  geom_sf(data = countries_sf %>% filter((name_en == "Indonesia")),alpha = 0, size=0.5,fill="white")+
  geom_sf(data = countries_sf %>% filter((name_en == "Malaysia") | (name_en == "Brunei")),size=0.5,alpha=0,fill="white") +
  coord_sf(xlim = c(109.25, 119.1), ylim = c(-4,7)) +
  scalebar(dist = 100, dist_unit = "km",location="bottomright",st.bottom=TRUE,height=0.75,st.size = 3.5,st.dist=0.75,border.size = 0.5,
  transform = TRUE, model = "WGS84",x.min=116.50,x.max=118.75,y.min=-4.15,y.max=-3.9)+
  geom_text(aes(x = 115.65, y = 1.1, label = "B",vjust = 0, hjust = -1),size=4) +
  geom_text(aes(x = 112.40, y =-0.85, label = "C",vjust = 0, hjust = -1),size=4) +
  geom_text(aes(x = 118.5, y = 6.92, label = "A",vjust = 0, hjust = -1),size=6) +
  scale_fill_brewer(palette = "RdBu", drop = FALSE, direction = -1,na.translate=F,
  aesthetics = c("colour", "fill"), name = "Change in deforestation\nprobability (pp)") +
  new_scale_fill()+
  geom_point(data=subset(mill_df,cert_status=="certified"),aes(x=longitude,y=latitude,fill=cert_status),size=1.2,shape=21,colour="grey20")+
  scale_fill_manual(name="",values="grey20",labels="Mill certified in 2016") +
  new_scale_fill() +
  geom_point(data=subset(mill_df,cert_status=="not-certified"),aes(x=longitude,y=latitude,fill=cert_status),size=1.2,shape=24,colour="grey20")+
  scale_fill_manual(name="",values="grey20",labels="Non certified Mill") +
  geom_rect(aes(xmin=115.45,xmax=116.45, ymin=-0.34,ymax=1.0), fill=NA,size=1.2,colour="black", alpha=0.2)+
  geom_rect(aes(xmin=112.25,xmax=113.25, ymin=-0.96,ymax=-2.30), fill=NA,size=1.2,colour="black", alpha=0.2)+
  theme(plot.margin=unit(c(1,0,0.5,0),"cm"),
  legend.margin = margin(-0.75,0,0,0, unit="cm")) +
  theme_map

oa_map

# create inset map B
inset_map1 <- ggplot(sample_df) +  
  geom_tile(aes(x=x, y=y, fill=value),alpha=0.8,size=0.0001,show.legend = FALSE) + 
  coord_sf() +
  geom_sf(data = countries_sf %>% filter((name_en == "Indonesia")),alpha = 0, size=0.5,fill="white")+
  geom_sf(data = countries_sf %>% filter((name_en == "Malaysia") | (name_en == "Brunei")),size=0.5,alpha=0,fill="white") +
  coord_sf(xlim = c(115.45, 116.45), ylim = c(-0.34,1.0)) +
  scalebar(dist = 25, dist_unit = "km",st.bottom=TRUE,height=0.2,st.size = 3,st.dist=0.25,border.size = 0.5,
  transform = TRUE, model = "WGS84",x.min=115.40,x.max=115.95,y.min=0.95,y.max=1.1)+
  scale_fill_brewer(palette = "RdBu", drop = FALSE, direction = -1,na.translate=F,
  aesthetics = c("colour", "fill"), name = "Change in\ndeforestation probability (%)") +
  new_scale_fill()+
  geom_point(data=subset(mill_df,cert_status=="certified"),aes(x=longitude,y=latitude,fill=cert_status),size=2,shape=21,colour="black",show.legend = FALSE)+
  scale_fill_manual(name="",values="grey20",labels="Mill certified in 2016") +
  new_scale_fill() +
  geom_point(data=subset(mill_df,cert_status=="not-certified"),aes(x=longitude,y=latitude,fill=cert_status),size=2,shape=24,colour="black",show.legend = FALSE)+
  scale_fill_manual(name="",values="grey20",labels="Non certified Mill") +
  geom_text(aes(x = 116.25, y = 0.92, label = "B",vjust = 0, hjust = -1),size=6) +
  geom_text(aes(x = 113.05, y = -1.05, label = "C",vjust = 0, hjust = -1),size=6) +
  theme_map +
  theme(plot.margin=unit(c(0,1.5,0.5,-7),"cm"))
 

inset_map1

# create inset map C
inset_map2 <- ggplot(sample_df) +  
  geom_tile(aes(x=x, y=y, fill=value,color=value),alpha=0.8,size=0.0001,show.legend = FALSE) + 
  coord_sf() +
  geom_sf(data = countries_sf %>% filter((name_en == "Indonesia")),alpha = 0, size=0.5,fill="white")+
  geom_sf(data = countries_sf %>% filter((name_en == "Malaysia") | (name_en == "Brunei")),size=0.5,alpha=0,fill="white") +
  coord_sf(xlim = c(112.25, 113.25), ylim = c(-0.96,-2.30)) +
  scalebar(dist = 25, dist_unit = "km",st.bottom=TRUE,height=0.075,st.size = 3,st.dist=0.1,border.size = 0.5,
  transform = TRUE, model = "WGS84",x.min=112.35,x.max=112.75,y.min=-0.60,y.max=-1.00)+
  scale_fill_brewer(palette = "RdBu", drop = FALSE, direction = -1,na.translate=F,
  aesthetics = c("colour", "fill"), name = "Change in\ndeforestation probability (%)") +
  new_scale_fill()+
  geom_point(data=subset(mill_df,cert_status=="certified"),aes(x=longitude,y=latitude,fill=cert_status),size=2,shape=21,colour="black",show.legend = FALSE)+
  scale_fill_manual(name="",values="grey20",labels="Mill certified in 2016") +
  new_scale_fill() +
  geom_point(data=subset(mill_df,cert_status=="not-certified"),aes(x=longitude,y=latitude,fill=cert_status),size=2,shape=24,colour="black",show.legend = FALSE)+
  scale_fill_manual(name="",values="grey20",labels="Non certified Mill") +
  guides(fill = FALSE, color = FALSE, linetype = FALSE, shape = FALSE) +
  geom_text(aes(x = 113.05, y = -1.05, label = "C",vjust = 0, hjust = -1),size=6) +
  #labs(title = "B") +
  theme_map +
  theme(plot.margin=unit(c(0,1.5,0.5,-7),"cm"))

inset_map2

# combine overview and inset maps
comb_defor_dif_map <- ggdraw(xlim = c(0, 35), ylim = c(0, 22)) +
  draw_plot(oa_map, x = 1, y = 0, width = 22, height = 22) +
  draw_plot(inset_map1, x = 29, y = 10, width = 10.5, height = 10.5) +
  draw_plot(inset_map2, x = 29, y = 0.0, width = 10.5, height = 10.5)

comb_defor_dif_map

# save to png file
ggsave(comb_defor_dif_map,file="output/fig05_comb_ioc_deforprob_map.png",dpi=500,w=8,h=6,unit="in",type="cairo-png")

#### Paper calculations ####
## Section 3.1
# We drew an evenly spaced, 2 x 2 km grid across the supply shed of all 
# Kalimantan oil palm mills and sampled individual points at the intersection 
# points of  this grid (106,889 points, representing 422,556 km2). 
n <- length(unique(full.df$sid))
area <- n * area_wght

## Section 3.4
# for the 11,261 points that fell within one of 890 oil palm concessions (17), 
# we assigned the oil palm concession as the management unit; (2) for the remaining 
# 22,301 points that fell within one of 396 timber and pulp concessions (46,47), we 
# assigned the timber or pulp concession as the management unit; (3) for the remaining 
# 25,951 points falling outside of any concession, we assigned the local village (Desa) 
op_npts <- df %>% 
  filter(year==2004, conc_class==3) %>% 
  n_distinct()
timber_npts <- df %>% 
  filter(year==2004, conc_class==2) %>% 
  n_distinct()
other_npts <- df %>% 
  filter(year==2004, conc_class==1) %>% 
  n_distinct()


## Section 4.1
# Total Kalimantan mill capacity increased 496% from 2,373 to 14,145 tonnes of FFB hour-1 
# between 2004 and 2016. Mill investments were particularly pronounced in Central and East 
# Kalimantan, where installed capacity increased 548% and 2,450%, respectively. The first 
# letter of intent to RSPO-certify a Kalimantan oil palm mill was issued on March 12, 2009. 
# Subsequently, certification rates increased in all provinces except North Kalimantan. By 
# 2016, 23% of total mill capacity in Kalimantan was certified, with an additional 25% 
#   attributed to mills held by RSPO member companies. 

# Get percent cap change from 2004 to 2016
total_mill_cap <- prov_caps_long_df %>%
  group_by(year) %>%
  summarise (totals=sum(value))

# Overall
pc_cap_change = (max(total_mill_cap$totals) -min(total_mill_cap$totals)) / min(total_mill_cap$totals) *100

# Get percent cap change from 2004 to 2016
total_mill_cap_prov <- prov_caps_long_df %>%
  group_by(prov,year) %>%
  summarise (totals=sum(value))

# Percentage by certification status (2016)
pc_cert <- prov_caps_long_df %>%
  group_by(cert) %>%
  filter(year == "2016-01-01") %>%
  summarise (totals=sum(value)) %>%
  mutate(pc_share = paste0(round(100 * totals/sum(totals), 0), "%"))

## Section 4.4
# Compared to the no-certification counterfactual, corporate group and local 
# supply chain spillovers generated a statistically insignificant increase in 
# total 2016 forest area outside of certified oil palm supply bases (median 
# simulation = +9 km2, 0.025-0.975 quantile range = -46 to +102 km2) (Figure 
# 4 and Figure 5). However, these effects differed across with land use zones. 
# Spillovers had a clear conservation benefit within the forest estate 
# (median simulation = +49 km2, 0.025-0.975 quantile range = +12 to +141 km2), 
# but induced deforestation in APL lands (median simulation = -44 km2, 0.025-0.975 
# quantile range = -73 to -15 km2). The median simulation indicates that aggregate 
# avoided deforestation from certification’s direct and spillover effects in 
# Kalimantan was 58 km2. However, due to the large variance in simulated spillover 
# effects, only 85% of simulations estimated a net increase in total Kalimantan 
# forest area. This contrasts with the precisely estimated, positive direct 
# impacts of certification in certified supply bases (median simulation = 
# +28 km2, 0.025-0.975 quantile range = +19 to +37 km2).
all.dif %>% 
  filter(op_type %in% c(1, 2), kh==-1) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE) %>% 
  quantile(c(0.025, 0.5, 0.975))

all.dif %>% 
  filter(op_type %in% c(1,2), kh==1) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE) %>% 
  quantile(c(0.025, 0.5, 0.975))

all.dif %>% 
  filter(op_type %in% c(1,2), kh==0) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE) %>% 
  quantile(c(0.025, 0.5, 0.975))

x <- all.dif %>% 
  filter(op_type == 0, kh==-1) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE)

all.dif %>% 
  filter(op_type == 0, kh==-1) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE) %>% 
  quantile(c(0.5))

1-ecdf(x)(0)

all.dif %>% 
  filter(op_type == 3, kh==-1) %>% 
  select(ends_with("cum_defor_defor")) %>% 
  unlist(., use.names = FALSE) %>% 
  quantile(c(0.025, 0.5, 0.975))

## 5.3
# We estimate that, since the first RSPO certificate was issued in 
# 2009, Kalimantan’s certified and non-certified oil palm concessions 
# experienced 33,000 km2 of forest loss. Our maximum estimate  of 
# avoided deforestation from the RSPO’s direct and indirect impacts 
# is <1% of this clearing. 
op_defor <- sample %>%
  filter(conc_class == 3,
         year >= 2009) %>%
  summarise(defor_area = sum(defor, na.rm = TRUE) * area_wght)

max_est <- all.dif %>% 
  filter(op_type == 0, kh==-1) %>% 
  max()

max_est / op_defor

## SI Section 9.2
# Land use zoning
short_df <- df %>%
  filter(year==2016) %>%
  mutate(kh_release = (kh_2005==1) & (kh_2018==0))

release_rate <- short_df %>% 
  filter(on_op==1, kh_2005==1) %>% 
  tabyl(kh_release, on_rspo)
(release_rate %>% filter(kh_release == 1)) / (release_rate %>% summarise_all(sum))











## DEPRECATED

# KH release and relation to deforestation
short_df <- df %>%
  filter(year==2016, op_type==1) %>%
  mutate(kh_release = (kh_2005==1) & (kh_2018==0),
         cert_shed = pccp>0)
release_rate <- short_df %>% tabyl(kh_release, cert_shed)
release_rate %>% filter(kh_release == TRUE) / release_rate %>% summarise_all(sum)

permitted_defor <- short_df %>%
  filter(kh05_defor==1) %>%
  tabyl(kh_release, on_rspo)
(permitted_defor %>% select(kh_release==TRUE))


sample <- sample %>%
  mutate(kh_release = (kh_2005==1) & (kh_2018==0),
         kh_all = kh_2005,
         kh_all = replace(kh_all, kh_release==1, 2))



## Test robustness to sample with only large management units
points_pid <- sample %>% 
  filter(year==2004) %>%
  group_by(pid) %>% 
  summarize(points_pid = n_distinct(sid))

sample <- sample %>%
  full_join(points_pid, by = "pid")

allkh_mod <- bife(defor ~ cert_now:factor(kh) + pccp:factor(kh) + g_cert_shr:factor(kh) + 
                    factor(kh) + factor(year) | pid, data = sample %>% filter(points_pid>10))
allkh_bc <- bias_corr(allkh_mod)
allkh_ape <- allkh_mod %>% bc_APE()



## Section 4.2 (Cut section)
cert_defor <- sample %>%
  filter(cert_now==1) %>%
  select(defor)
defor <- cert_defor %>%
  filter(defor==TRUE) %>%
  dim()
nodefor <- cert_defor %>%
  filter(defor==FALSE) %>%
  dim()
defor_rate <- defor[1] / (defor[1]+nodefor[1])
att <- -0.0119
cf <- defor_rate - att
pct_change <- att / cf
