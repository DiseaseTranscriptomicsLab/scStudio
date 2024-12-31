###############################################
#
# UI for homepage tab
#
###############################################

tab_HOME <- tabPanel(
  title = "",
  icon = NULL,
    #div(class="fa fa-home", role = "navigation"), 
  value = "home", # tab ID
 
  htmlTemplate("www/landing-page.html",
               latest_updates_button = actionButton('btn_indicator_updates', 
                                                    "View recent updates", class = "button")
            
  ))



