# S1 still needs an app mount

Neither `IE_CROSSMODAL_CARDS` in `apps/islet/R/mod_crossmodal_ui.R` nor
`IE_LYMPHNODE_CARDS` in `apps/islet/R/mod_lymphnode_ui.R` lists
`s1_size_composition`. No other app module calls its UI or server, so the current
app cannot display S1 even after correcting the card's active-tab guard.

The card now defaults to `tab_key = "Plot"`, the valid `IE_PIPELINE` key whose
subtitle is "Islet size · composition". This is a proposed destination, not an
existing mount. The integration owner needs to mount
`card_s1_size_composition_ui(ns)` and `card_s1_size_composition_server(...)`
together on that tab, forwarding `active_tab = reactive(input$tabs)`. If the owner
chooses another tab, pass its exact key as `tab_key`.

Registry and app wiring changes are outside this worker's file ownership in
`BRIEF.md`; `FIX.md` explicitly requests this handoff when S1 is unregistered.
