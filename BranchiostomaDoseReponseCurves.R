#Author: Luis A. Yanez-Guerra 
# load packages -----------------------------------------------------------
library(pacman)
p_load(tidyverse, png, patchwork, drc, devtools, vroom, credentials, ggsci)

#check working directory- should be the project dir, all directories will be defined relative to this dir
Sys.setlocale("LC_MESSAGES", 'en_GB.UTF-8')
Sys.setenv(LANG = "en_US.UTF-8")
set.seed(666)
# read data ---------------------------------------------------------------
#Allmixturesofpeptides
curve <- vroom("SupplementaryFile11_rawdata_receptor_characterization.csv")

#convert to tibble Dose-Responsecurves
AllGPCRtoplot <- curve %>%
  pivot_longer(starts_with(c("0", "1e", "1", "Neg")), 
               names_to = "concentration", values_to = "luminescence")


#convert conc valus to double
AllGPCRtoplot$concentration <- as.double(AllGPCRtoplot$concentration)

# Delete the N.A ----------------------------------------------------------
AllGPCRtoplot <- AllGPCRtoplot[!is.na(AllGPCRtoplot$luminescence), ]                 # Omit NA by column via is.na

# Functions ---------------------------------------------------------------
#FUNCTION to normalise to reference (zero ligand control)
normalize_to_ctr <- function(x) {
  return (100*(x - x[1]) / (max(x) - x[1]))
}



#Function to add grids to the curves (THIS CODE IS FROM https://github.com/kassambara/ggpubr/blob/master/R/grids.R)
#Thanks to kassambara for making it public and share it!.

grids <- function(axis = c("xy", "x", "y"), color = "slategray1", size = NULL, linetype = NULL)
{
  axis <- match.arg(axis)
  grid.major <- element_line(color = color, size = size,
                             linetype = linetype)
  grid.minor <- element_line(color = color, size = 0.25,
                             linetype = linetype)
  
  switch(axis,
         xy = theme(panel.grid.major = grid.major, panel.grid.minor = grid.minor),
         x = theme(panel.grid.major.x = grid.major, panel.grid.minor.x = grid.minor),
         y = theme(panel.grid.major.y = grid.major, panel.grid.minor.y = grid.minor)
  )
}


#This is normalization of the receptors by the plate.
AllGPCRtoplotTR <- AllGPCRtoplot %>%
  group_by(Transfection, Receptor)%>%
  mutate('norm_luminescence'=normalize_to_ctr(luminescence))


# Plot with boxplot -------------------------------------------------------
#This helps to plot the dose-response curve normalized.

AllGPCRtoplotTR %>% 
  ggplot(aes(x=concentration, y=norm_luminescence, colour=Peptide, group=Peptide)) +
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE) +
  scale_x_log10(breaks = c(1e-13, 1e-11,1e-9,1e-7,1e-5, 1e-3),limits = c(1e-12,1e-3)) +
  theme_half_open() +
  scale_color_aaas() +
  theme(axis.text=element_text(size = 12), 
        legend.text = element_text(size=12), 
        legend.title=element_text(size=12),
        axis.title=element_text(size=20), 
        axis.title.x=element_text(margin = margin(t = 19), (face = "bold")),
        panel.background = element_blank(),
        legend.position = "bottom") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4) +
  stat_summary(fun.data = mean_se, geom="errorbar", size=0.3) +
  facet_wrap(vars(Receptor)) + 
  grids(linetype = "longdash")  #to change the colour of grids, change it in the function above 


#This helps to plot the raw data for comparison, notice that the Corazonin receptors produce
#a much lower luminescence after the addition of peptides than the GnRH receptors.
#However, once normalized, you can see the curve

AllGPCRtoplotTR %>% 
  ggplot(aes(x=concentration, y=luminescence, colour=Peptide, group=Peptide)) +
  geom_smooth(method = drm, method.args = list(fct = L.4()), se = FALSE) +
  scale_x_log10(breaks = c(1e-11,1e-9,1e-7,1e-5, 1e-3),limits = c(1e-12,1e-3)) +
  theme_half_open() +
  scale_color_aaas() +
  theme(axis.text=element_text(size = 12), 
        legend.text = element_text(size=12), 
        legend.title=element_text(size=12),
        axis.title=element_text(size=20), 
        axis.title.x=element_text(margin = margin(t = 19), (face = "bold")),
        panel.background = element_blank(),
        legend.position = "bottom") +
  stat_summary(fun.y=mean, geom="point", shape=20, size=4) +
  stat_summary(fun.data = mean_se, geom="errorbar", size=0.3) +
  facet_wrap(vars(Receptor)) + 
  grids(linetype = "longdash")  #to change the colour if grids, change it in the function above 

