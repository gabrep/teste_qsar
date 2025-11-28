library(tidyverse)
setwd('r/')
pred.test <- read.csv('../data/Prediction_models_test.csv')
pred.train <- read.csv('../data/Prediction_models_train.csv')

pred.train %>% reorder(Model, R.Squared, decreasing = F)

pred.train %>% 
  arrange(desc(R.Squared)) %>% 
  mutate(grupo_fill = ifelse(Model == "RandomForestRegressor",
                             "RF", "Outros")) %>% 
  slice(1:10) %>% 
  ggplot(aes(R.Squared, y = reorder(Model, R.Squared,decreasing = F), fill=grupo_fill))+
  geom_col()+
  scale_fill_manual(values=c('RF' = 'deepskyblue', 'Outros' = 'grey90'), guide='none')+
  labs(y=NULL)+
  theme_classic()+
  geom_text(aes(label = round(R.Squared, 3)), hjust=1.1, size=3)+
  coord_cartesian(xlim = c(0,1.1), expand=F)+
  theme(axis.text = element_text(color='black'))
ggsave('../figures/Pred_train_R.png', width = 4, height = 3)


pred.test %>% 
  arrange(desc(R.Squared)) %>% 
  mutate(grupo_fill = ifelse(Model == "RandomForestRegressor",
                             "RF", "Outros")) %>% 
  slice(1:10) %>% 
  ggplot(aes(R.Squared, y = reorder(Model, R.Squared,decreasing = F), fill=grupo_fill))+
  geom_col()+
  scale_fill_manual(values=c('RF' = 'deepskyblue', 'Outros' = 'grey90'), guide='none')+
  labs(y=NULL)+
  theme_classic()+
  geom_text(aes(label = round(R.Squared, 3)), hjust=1.1, size=3)+
  coord_cartesian(xlim = c(0,1), expand=F)+
  theme(axis.text = element_text(color='black'))
ggsave('../figures/Pred_test_R.png', width = 4, height = 3)

lotus <- read.csv('../data/Lotus_predito_top10.csv')

lotus %>%
  mutate(IC50_nM = IC50 * 1000000000) %>% 
  ggplot(aes(IC50_nM, reorder(X1, IC50, decreasing = F), fill=IC50_nM))+
  geom_col()+
  geom_text(aes(label = paste('MW = ', MW, 'LogP = ', round(LogP,2))), size=3, hjust=c(rep(-.01, 3), rep(1.05,7)))+
  labs(y=NULL)+
  scale_fill_gradient(low = 'pink' , high = 'violetred')+
  theme_classic()+
  coord_cartesian(xlim=c(0,155), expand = F)+
  theme(axis.text = element_text(color='black'))
ggsave('../figures/Lotus_top10.png', width = 6, height = 5)
