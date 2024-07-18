
fit_geohisse= function(phy, data, trans_simple, trans_hidden, sann = F){
  ### list to keep all models
  model_list = list()
  ### model 0 - dispersal only
  model_list$model_0 = GeoHiSSE(
    phy = phy, 
    data = data, 
    f= sampling_f, 
    turnover=c(1,1,0), 
    eps=c(1,1), 
    hidden.states= FALSE,
    trans.rate = trans_simple, 
    root.type="madfitz",
    sann=sann
  )
  ### model 1 - range dependent
  model_list$model_1 = GeoHiSSE(
    phy = phy, 
    data =data, 
    f= sampling_f, 
    turnover=c(1,2,3), 
    eps=c(1,1), 
    hidden.states= FALSE,
    trans.rate = trans_simple, 
    root.type="madfitz",
    sann=sann
  )
  ### model 2 - hidden diversification
  model_list$model_2 = GeoHiSSE(
    phy = phy, 
    data = data, 
    f= sampling_f, 
    turnover=c(1,1,0,2,3,0), 
    eps=c(1,1,1,1), 
    hidden.states= T,
    trans.rate = trans_hidden, 
    root.type="madfitz",
    sann=sann
  )
  ### model 3 - range dependent with hidden state
  model_list$model_3 = GeoHiSSE(
    phy = phy, 
    data = data, 
    f= sampling_f, 
    turnover=c(1,2,3,4,5,6), 
    eps=c(1,1,1,1), 
    hidden.states= T,
    trans.rate = trans_hidden, 
    root.type="madfitz",
    sann=sann
  )
  ### AICc weights
  aicc_w = GetAICWeights(list(model_0 = model_list$model_0,
                              model_1 = model_list$model_1,
                              model_2 = model_list$model_2,
                              model_3 = model_list$model_3
  ),
  criterion="AICc")
  ### name of the best model 
  best_model_name = names(aicc_w[aicc_w == max(aicc_w)])
  ### pick the best model
  best_model = model_list[best_model_name]
  ### return
  return(best_model)
}
