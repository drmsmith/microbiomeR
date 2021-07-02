library(ggplot2)
library(cowplot)
library(ggpubr)

source('solve.R')

##################################
### FIGURE 1: Model Comparison ###
##################################

cols_interactions = c("#1b9e77","#7570b3","#e7298a")
labs_interactions = c(expression(paste(epsilon)), expression(paste(eta)), expression(paste(phi)))

### Default parameters (r_R = 0.8)

# coordinates for max prevalence
max_C_R_x_model1a = 0
max_C_R_y_model1a = fig1_model1a$prevalence_C_R[which(fig1_model1a$prevalence_C_R == max(fig1_model1a$prevalence_C_R))]

max_C_R_x_model2a = fig1_model2a$par_val[which(fig1_model2a$prevalence_C_R == max(fig1_model2a$prevalence_C_R))]
max_C_R_y_model2a = fig1_model2a$prevalence_C_R[which(fig1_model2a$prevalence_C_R == max(fig1_model2a$prevalence_C_R))]

max_C_R_x_model3a_epsilon = fig1_model3a_epsilon_med$par_val[which(fig1_model3a_epsilon_med$prevalence_C_R == max(fig1_model3a_epsilon_med$prevalence_C_R))]
max_C_R_y_model3a_epsilon = fig1_model3a_epsilon_med$prevalence_C_R[which(fig1_model3a_epsilon_med$prevalence_C_R == max(fig1_model3a_epsilon_med$prevalence_C_R))]

max_C_R_x_model3a_eta = fig1_model3a_eta_med$par_val[which(fig1_model3a_eta_med$prevalence_C_R == max(fig1_model3a_eta_med$prevalence_C_R))]
max_C_R_y_model3a_eta = fig1_model3a_eta_med$prevalence_C_R[which(fig1_model3a_eta_med$prevalence_C_R == max(fig1_model3a_eta_med$prevalence_C_R))]

max_C_R_x_model3a_phi = fig1_model3a_phi_med$par_val[which(fig1_model3a_phi_med$prevalence_C_R == max(fig1_model3a_phi_med$prevalence_C_R))]
max_C_R_y_model3a_phi = fig1_model3a_phi_med$prevalence_C_R[which(fig1_model3a_phi_med$prevalence_C_R == max(fig1_model3a_phi_med$prevalence_C_R))]

max_C_R_x_model3a = c(max_C_R_x_model3a_epsilon, max_C_R_x_model3a_eta, max_C_R_x_model3a_phi)
max_C_R_y_model3a = c(max_C_R_y_model3a_epsilon, max_C_R_y_model3a_eta, max_C_R_y_model3a_phi)

ymax_a = max(c(max_C_R_y_model1a, max_C_R_y_model2a, max(fig1_model3a_interactions$high))) + 0.02

# plot model 1
p_model1a = ggplot(fig1_model1a, aes(x = par_val, y = prevalence_C_R, colour = time))+
  geom_line(stat = 'identity', colour = 'black')+
  ylab('Pathogen colonization prevalence')+xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme_classic()+
  # arrow lines
  geom_segment(x = max_C_R_x_model1a, xend = max_C_R_x_model1a, y = max_C_R_y_model1a, yend = -0.015,
               linetype = 2,colour = 'black', size = 0.2, show.legend = F)+
  # arrow heads
  geom_segment(x = max_C_R_x_model1a, xend = max_C_R_x_model1a, y = -0.014, yend = -0.015,
               linetype = 1, colour = 'black', size = 0.2, arrow = arrow(length = unit(0.1, "inches")), show.legend = F)+
  scale_colour_manual(values = c('black'), name = 'Strain', labels = expression(paste(C[r])))+
  theme(legend.position = 'bottom')+
  ylim(0, ymax_a)
p_model1a

# plot model 2
p_model2a = ggplot(fig1_model2a_strains, aes(x = par_val, y = value, colour = strain, alpha = strain, linetype = strain))+
  # arrow lines
  geom_segment(x = max_C_R_x_model2a, xend = max_C_R_x_model2a, y = max_C_R_y_model2a, yend = -0.015,
               linetype = 2,colour = 'red', size = 0.2, show.legend = F)+
  # arrow heads
  geom_segment(x = max_C_R_x_model2a, xend = max_C_R_x_model2a, y = -0.014, yend = -0.015,
               linetype = 1, colour = 'red', size = 0.2, arrow = arrow(length = unit(0.1, "inches")), show.legend = F)+
  geom_line(stat = 'identity')+
  scale_colour_manual(values = c('red', 'black'), name = 'Strain', labels = c('Resistant', 'Sensitive'))+
  scale_alpha_manual(values = c(1,0.5), name = 'Strain', labels = c('Resistant', 'Sensitive'))+
  scale_linetype_manual(values = c(1,1), name = 'Strain', labels = c('Resistant', 'Sensitive'))+
  ylab('Pathogen colonization prevalence')+xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme_classic()+
  theme(legend.position = 'bottom')+
  ylim(0, ymax_a)
p_model2a

# plot model 3
p_model3a = ggplot(fig1_model3a_interactions, aes(x = par_val, y = medium, colour = interaction))+
  scale_colour_manual(values = cols_interactions, name = 'Interaction', labels = labs_interactions)+
  # arrow lines
  geom_segment(x = max_C_R_x_model3a[1], xend = max_C_R_x_model3a[1], y = max_C_R_y_model3a[1], yend = -0.015,
               linetype = 2, colour = cols_interactions[1], size = 0.2)+
  geom_segment(x = max_C_R_x_model3a[2], xend = max_C_R_x_model3a[2], y = max_C_R_y_model3a[2], yend = -0.015,
               linetype = 2, colour = cols_interactions[2], size = 0.2)+
  geom_segment(x = max_C_R_x_model3a[3], xend = max_C_R_x_model3a[3], y = max_C_R_y_model3a[3], yend = -0.015,
               linetype = 2, colour = cols_interactions[3], size = 0.2)+
  # # arrow heads
  geom_segment(x = max_C_R_x_model3a[1], xend = max_C_R_x_model3a[1], y = -0.014, yend = -0.015,
               linetype = 1, colour = cols_interactions[1], size = 0.2, arrow = arrow(length = unit(0.1, "inches")))+
  geom_segment(x = max_C_R_x_model3a[2], xend = max_C_R_x_model3a[2], y = -0.014, yend = -0.015,
               linetype = 1, colour = cols_interactions[2], size = 0.2, arrow = arrow(length = unit(0.1, "inches")))+
  geom_segment(x = max_C_R_x_model3a[3], xend = max_C_R_x_model3a[3], y = -0.014, yend = -0.015,
               linetype = 1, colour = cols_interactions[3], size = 0.2, arrow = arrow(length = unit(0.1, "inches")))+
  geom_line(stat = 'identity')+
  theme_classic()+
  ylab('Pathogen colonization prevalence')+xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme_classic()+
  geom_ribbon(aes(ymin = low, ymax = high, fill = interaction), alpha = 0.1, colour = NA)+
  scale_fill_manual(values = cols_interactions, name = 'Interaction', labels = labs_interactions)+
  theme(legend.position = 'bottom')+
  ylim(0, ymax_a)
p_model3a

### Perfect resistance (r_R = 1)

ymax_b = max(c(fig1_model1b$prevalence_C_R, fig1_model2b_strains$value, fig1_model3b_interactions$high)) + 0.02

# plot model 1
p_model1b = ggplot(fig1_model1b, aes(x = par_val, y = prevalence_C_R, colour = time))+
  geom_line(stat = 'identity', colour = 'black')+
  ylab('Pathogen colonization prevalence')+xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme_classic()+
  scale_colour_manual(values = c('black'), name = 'Strain', labels = expression(paste(C[r])))+
  theme(legend.position = 'bottom')+
  ylim(0, ymax_b)
p_model1b

# plot model 2
p_model2b = ggplot(fig1_model2b_strains, aes(x = par_val, y = value, colour = strain, alpha = strain, linetype = strain))+
  geom_line(stat = 'identity')+
  scale_colour_manual(values = c('red', 'black'), name = 'Strain', labels = c('Resistant', 'Sensitive'))+
  scale_alpha_manual(values = c(1,0.5), name = 'Strain', labels = c('Resistant', 'Sensitive'))+
  scale_linetype_manual(values = c(1,1), name = 'Strain', labels = c('Resistant', 'Sensitive'))+
  ylab('Pathogen colonization prevalence')+xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme_classic()+
  theme(legend.position = 'bottom')+
  ylim(0, ymax_b)
p_model2b

# plot model 3
p_model3b = ggplot(fig1_model3b_interactions, aes(x = par_val, y = medium, colour = interaction))+
  scale_colour_manual(values = cols_interactions, name = 'Interaction', labels = labs_interactions)+
  geom_line(stat = 'identity')+
  theme_classic()+
  ylab('Pathogen colonization prevalence')+xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme_classic()+
  geom_ribbon(aes(ymin = low, ymax = high, fill = interaction), alpha = 0.1, colour = NA)+
  scale_fill_manual(values = cols_interactions, name = 'Interaction', labels = labs_interactions)+
  theme(legend.position = 'bottom')+
  ylim(0, ymax_b)
p_model3b


# combine plots and save
p_fig1 = plot_grid(p_model1a, p_model2a+ylab(' '), p_model3a+ylab(' '), ncol = 3, nrow = 1, align = 'h', axis = 'b')
p_fig1
ggsave(p_fig1, width = 10, height = 4.2, filename = paste0(getwd(),"/figures/Fig1.png"))
ggsave(p_fig1, width = 10, height = 4.2, filename = paste0(getwd(),"/figures/Fig1.pdf"))

p_fig1_perfectR = plot_grid(p_model1b, p_model2b+ylab(' '), p_model3b+ylab(' '), ncol = 3, nrow = 1, align = 'h', axis = 'b')
p_fig1_perfectR
ggsave(p_fig1_perfectR, width = 10, height = 4.2, filename = paste0(getwd(),"/figures/Fig1_perfectR.png"))
ggsave(p_fig1_perfectR, width = 10, height = 4.2, filename = paste0(getwd(),"/figures/Fig1_perfectR.pdf"))

################
### FIGURE 2 ###
################

### Default parameters
p_fig2_final = ggplot(fig2_default, aes(x = par1_val, y = par2_val, colour = R_rate, size = prevalence_C_R))+
  geom_point(stat='identity')+
  theme_classic()+
  scale_colour_distiller(palette = 'Spectral', name = 'Resistance rate', limits = c(min(fig2_default$R_rate), max(fig2_default$R_rate)))+
  scale_size_continuous(name = expression(paste('Pathogen prevalence')), breaks = seq(0.11,0.6,0.02), range = c(0,10))+
  theme(text = element_text(size = 14),
        axis.text = element_text(size = 10))+
  xlab(expression(paste('Rate of antibiotic-induced pathogen clearance  ', (theta[c]))))+
  ylab(expression(paste('Rate of antibiotic-induced microbiome dysbiosis  ', (theta[m]))))
p_fig2_final

### Grid over range of parameters
p_fig2_mixed = ggplot(fig2_mixed, aes(x = par1_val, y = par2_val, colour = R_rate, size = prevalence_C_R))+
  geom_point(stat='identity')+
  theme_classic()+
  scale_colour_distiller(palette = 'Spectral', name = 'Resistance rate', limits = c(min(fig2_mixed$R_rate), max(fig2_mixed$R_rate)))+
  scale_size_continuous(name = 'Pathogen prevalence', limits = c(min(fig2_mixed$prevalence_C_R), max(fig2_mixed$prevalence_C_R)), range = c(-1,4), breaks = seq(0.05,0.8,0.05))+
  theme(text = element_text(size = 16),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 10))+
  xlab(expression(paste('Rate of antibiotic-induced pathogen clearance  ', (theta[c]))))+
  ylab(expression(paste('Rate of antibiotic-induced microbiome dysbiosis  ', (theta[m]))))+
  facet_wrap(facets = vars(label))
p_fig2_mixed

## save plots
ggsave(p_fig2_final, width = 9.6, height = 7.5, filename = paste0(getwd(),'/figures/Fig2_final.pdf'))
ggsave(p_fig2_final, width = 9.6, height = 7.5, filename = paste0(getwd(),'/figures/Fig2_final.png'))

ggsave(p_fig2_mixed, width = 13, height = 13, filename = paste0(getwd(),'/figures/Fig2_grid.pdf'))
ggsave(p_fig2_mixed, width = 13, height = 13, filename = paste0(getwd(),'/figures/Fig2_grid.png'))

################
### FIGURE 3 ###
################

label_a = expression(paste('Antibiotic exposure prevalence (',italic(a),')'))
label_prevalence = expression(paste('Colonization prevalence  ',C^R))
label_incidence = expression(paste('Daily colonization incidence (% of patients)'))
label_R_rate = expression(paste('Resistance rate  ',C^R,'/(',C^S + C^R,')'))

### Default pars
fig3_prev = ggplot(fig3, aes(x = par_val, y = prevalence_C_R, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_prevalence)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_prev

fig3_inc = ggplot(fig3, aes(x = par_val, y = incidence_C_R_daily*100, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_incidence)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_inc

fig3_R_rate = ggplot(fig3, aes(x = par_val, y = R_rate, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_R_rate)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_R_rate


fig3_prev_R_rate = ggarrange(fig3_prev+xlab(''), fig3_R_rate, nrow = 2, common.legend = T, legend = 'right', align = 'v')
fig3_prev_R_rate


### Medium r_R
fig3_prev_medR = ggplot(fig3_medR, aes(x = par_val, y = prevalence_C_R, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_prevalence)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_prev_medR

fig3_inc_medR = ggplot(fig3_medR, aes(x = par_val, y = incidence_C_R_daily*100, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_incidence)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_inc

fig3_R_rate_medR = ggplot(fig3_medR, aes(x = par_val, y = R_rate, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_R_rate)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_R_rate_medR


### Low r_R
fig3_prev_lowR = ggplot(fig3_lowR, aes(x = par_val, y = prevalence_C_R, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_prevalence)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_prev_lowR

fig3_inc_lowR = ggplot(fig3_lowR, aes(x = par_val, y = incidence_C_R_daily*100, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_incidence)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_inc

fig3_R_rate_lowR = ggplot(fig3_lowR, aes(x = par_val, y = R_rate, colour = Interactions, linetype = HGT)) + 
  geom_line()+
  ylab(label_R_rate)+
  xlab(label_a)+
  scale_colour_manual(values = c('#ff7f00','#6a3d9a'))+
  #scale_alpha_manual(values = c(0.3,0.6,1))+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')))+
  theme_classic()+
  guides(colour = F)+
  theme(text = element_text(size = 12),
        axis.text = element_text(size = 10))
fig3_R_rate_lowR


fig3_prev_R_rate_lowR = ggarrange(fig3_prev_lowR+xlab(''), 
                                  fig3_R_rate_lowR, 
                                  nrow = 2, common.legend = T, legend = 'right', align = 'v')
fig3_prev_R_rate_lowR

### Show across 3 levels of resistance

fig3_prev_R_rate_acrossR = ggarrange(fig3_prev_lowR+
                                       labs(subtitle=expression(paste('Low resistance level (',r[R],'=0.2)')),
                                            x = element_blank())+
                                       theme(plot.subtitle = element_text(hjust = 0.5, face = 'bold'))+
                                       ylim(0.03,0.42)+
                                       geom_hline(yintercept = 0.05, colour = 'grey', alpha = 0.5, size = 2),
                                     fig3_prev_medR+
                                       labs(subtitle=expression(paste('Intermediate resistance level (',r[R],'=0.5)')),
                                            x = element_blank(), y = element_blank())+
                                       theme(plot.subtitle = element_text(hjust = 0.5, face = 'bold'))+
                                       ylim(0.03,0.42)+
                                       geom_hline(yintercept = 0.05, colour = 'grey', alpha = 0.5, size = 2), 
                                     fig3_prev+
                                       labs(subtitle=expression(paste('High resistance level (',r[R],'=0.8)')),
                                            x = element_blank(), y = element_blank())+
                                       theme(plot.subtitle = element_text(hjust = 0.5, face = 'bold'))+
                                       ylim(0.03,0.42)+
                                       geom_hline(yintercept = 0.05, colour = 'grey', alpha = 0.5, size = 2),
                                     fig3_R_rate_lowR+
                                       ylim(0.25,1)+
                                       geom_hline(yintercept = 0.5, colour = 'grey', alpha = 0.5, size = 2),
                                     fig3_R_rate_medR+
                                       labs(y = element_blank())+
                                       ylim(0.25,1)+
                                       geom_hline(yintercept = 0.5, colour = 'grey', alpha = 0.5, size = 2), 
                                     fig3_R_rate+
                                         labs(y = element_blank())+
                                         ylim(0.25,1)+
                                       geom_hline(yintercept = 0.5, colour = 'grey', alpha = 0.5, size = 2),
                                     ncol = 3, nrow = 2, common.legend = T, legend = 'right', align = 'hv')
fig3_prev_R_rate_acrossR

ggsave(fig3_prev_R_rate_acrossR, width = 12, height = 7, filename = paste0(getwd(),'/figures/Fig3_prev_R_rate_acrossR.pdf'))
ggsave(fig3_prev_R_rate_acrossR, width = 12, height = 7, filename = paste0(getwd(),'/figures/Fig3_prev_R_rate_acrossR.png'))

fig3_inc_acrossR = ggarrange(fig3_inc_lowR+
                               labs(subtitle=expression(paste('Low resistance level (',r[R],'=0.2)')))+
                               theme(plot.subtitle = element_text(hjust = 0.5, face = 'bold'))+
                               ylim(1,6.5),
                             fig3_inc_medR+
                               labs(subtitle=expression(paste('Intermediate resistance level (',r[R],'=0.5)')))+
                               theme(plot.subtitle = element_text(hjust = 0.5, face = 'bold'))+
                               labs(y = element_blank())+
                               ylim(1,6.5), 
                             fig3_inc+
                               labs(subtitle=expression(paste('High resistance level (',r[R],'=0.8)')))+
                               theme(plot.subtitle = element_text(hjust = 0.5, face = 'bold'))+
                               labs(y = element_blank())+
                               ylim(1,6.5),
                             ncol = 3,common.legend = T, legend = 'right', align = 'hv')
                                     
fig3_inc_acrossR
ggsave(fig3_inc_acrossR, width = 12, height = 4, filename = paste0(getwd(),'/figures/Fig3_inc_acrossR.pdf'))
ggsave(fig3_inc_acrossR, width = 12, height = 4, filename = paste0(getwd(),'/figures/Fig3_inc_acrossR.png'))

#################
### FIGURE S5 ###
#################

figDynamics_a = ggplot(dynamics_combined, aes(x = time, y = prevalence_C_R, colour = Interaction))+
  geom_line()+
  theme_classic()+
  ylab(expression(paste('Colonization prevalence  ',C^R)))+
  xlab('Time (days)')+
  geom_vline(xintercept = c(90,180, 270), size = 2, colour = 'darkgrey', alpha = 0.6)+
  scale_colour_manual("Microbiome interaction", values = c("#666666", cols_interactions, '#d95f02', '#a6761d'))+
  annotate('text', x = 95, y = 0.27, label = 'intervention 1:\nadjust\nprescribing 1', hjust = 0, size = 2.5)+
  annotate('text', x = 185, y = 0.27, label = 'intervention 2:\nadjust\nprescribing 2', hjust = 0, size = 2.5)+
  annotate('text', x = 275, y = 0.27, label = 'intervention 3:\nreduce\nprescribing', hjust = 0, size = 2.5)+
  ylim(min(dynamics_combined$prevalence_C_R), max(dynamics_combined$prevalence_C_R)+0.02)
figDynamics_a

figDynamics_b = ggplot(dynamics_combined, aes(x = time, y = R_rate, colour = Interaction))+
  geom_line()+
  theme_classic()+
  ylab(expression(paste('Resistance rate  ',C^R,'/(',C^S + C^R,')')))+
  xlab('Time (days)')+
  geom_vline(xintercept = c(90,180, 270), size = 2, colour = 'darkgrey', alpha = 0.6)+
  scale_colour_manual("Microbiome interaction", values = c("#666666", cols_interactions, '#d95f02', '#a6761d'))+
  annotate('text', x = 95, y = 0.6, label = 'intervention 1:\nadjust\nprescribing 1', hjust = 0, size = 2.5)+
  annotate('text', x = 185, y = 0.6, label = 'intervention 2:\nadjust\nprescribing 2', hjust = 0, size = 2.5)+
  annotate('text', x = 275, y = 0.6, label = 'intervention 3:\nreduce\nprescribing', hjust = 0, size = 2.5)+
  ylim(min(dynamics_combined$R_rate), max(dynamics_combined$R_rate)+0.04)
figDynamics_b

figDynamics = ggarrange(figDynamics_a+xlab(''), figDynamics_b, align = 'v', nrow = 2, common.legend = T, legend = 'right')
figDynamics

ggsave(figDynamics, width = 8, height = 6, filename = paste0(getwd(),'/figures/FigDynamics.png'))
ggsave(figDynamics, width = 8, height = 6, filename = paste0(getwd(),'/figures/FigDynamics.pdf'))


#################
### FIGURE S7 ###
#################

# a: compare parameters
p_hgt_compare_pars = figS7_a_final%>%
  dplyr::select(-c(none, low, high))%>%
  pivot_longer(-c(par, par_val), values_to = "difference", names_to = "HGT_rate")%>%
  ggplot(aes(x = par_val, y = difference, linetype = factor(HGT_rate)))+
  geom_line(stat= 'identity')+
  facet_wrap(facets = vars(par), scales = 'free_x', ncol = 3)+
  scale_linetype_manual(values = c(3:1), name = expression(paste('HGT rate (',chi,')')), labels = c('high', 'low', 'null'))+
  theme_classic()+
  xlab('Parameter value')+
  ylab(expression(atop(paste('Colonization prevalence (',C^R,'),'),'absolute difference from HGT-free baseline')))+
  theme(axis.text = element_text(size = 6))
p_hgt_compare_pars

# b: vary chi_d relative to chi_e
p_hgt_varyHGT = figS7_varyHGT%>%
  ggplot(aes(x = par_val, y = prevalence_C_R, colour = factor(ratio)))+
  geom_line(stat= 'identity')+
  scale_colour_manual(name = expression(paste(chi[d],'/',chi[e])), values = c('#dadaeb','#bcbddc','#9e9ac8','#756bb1','#54278f'))+
  theme_classic()+
  xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  ylab(expression(paste('Colonization prevalence (',C^R,')')))+
  theme(axis.text = element_text(size = 6))
p_hgt_varyHGT

# c: 
p_hgt_vary_c = figS7_c_vary%>%
  ggplot(aes(x = par_val, y = prev_difference, colour = factor(cost)))+
  geom_line(stat= 'identity')+
  scale_colour_manual(name = 'fitness cost\nof resistance', values = rev(c('#d73027','#fc8d59','#fee090','#91bfdb','#4575b4')))+
  theme_classic()+
  ylab(expression(atop(paste('Total colonization prevalence (',C^S + C^R,'),'),'absolute difference from HGT-free baseline')))+
  xlab(expression(paste('Antibiotic exposure prevalence (',italic(a),')')))+
  theme(axis.text = element_text(size = 6))+
  geom_hline(yintercept = 0, linetype = 2)
p_hgt_vary_c

p_hgt_sensitivity_grid = plot_grid(p_hgt_compare_pars, 
                                   plot_grid(p_hgt_varyHGT, p_hgt_vary_c, nrow = 2, align = 'v', axis = 'lr', labels = c('B', 'C')),
                                   ncol = 2, rel_widths = c(1.35,1), labels = c('A',''))
p_hgt_sensitivity_grid

ggsave(p_hgt_sensitivity_grid, width = 12.5, height = 8, filename = paste0(getwd(),'/figures/FigHGT.png'))
ggsave(p_hgt_sensitivity_grid, width = 12.5, height = 8, filename = paste0(getwd(),'/figures/FigHGT.pdf'))



