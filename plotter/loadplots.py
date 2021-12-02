import numpy as np
from html import add_web_section, PlotsInPipeline

def loadPlots(web, siminfo):

    PlotsInWeb = PlotsInPipeline()

    # Adding plots
    title = "Cosmic Scatter Rate"
    caption = ""
    filename = "Cosmic_scatter_rate.png"
    id = abs(hash("ScatterRate"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    # Setting section
    title = "Cosmic Scatter Rate"
    id = abs(hash("CosmicScatterRate"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    # Adding plots
    title = "Halo Mass Function"
    caption = ""
    filename = "HMF_%04d.png" % siminfo.n_snapshots
    id = abs(hash("HMF"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Subhalo Mass Function"
    caption = ""
    filename = "SubHMF_%04d.png" % siminfo.n_snapshots
    id = abs(hash("SubHMF"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    # Setting section
    title = "Halo/Subhalo Mass Functions"
    id = abs(hash("Halo/SubhaloMF"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    # Adding plots
    title = "Concentration-Mass Relation"
    caption = ""
    filename = "cM_relation.png"
    id = abs(hash("CM relation"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Vmax-Mass Relation"
    caption = ""
    filename = "VmaxM_relation.png"
    id = abs(hash("Vmax relation"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    # Setting section
    title = "Concentration/Vmax relations"
    id = abs(hash("cVmax"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    # Adding plots
    title = "Cross section-velocity relation"
    caption = ""
    filename = "particle_data_1.png"
    id = abs(hash("Particle data 1"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    # Setting section
    title = "Cross-section relations"
    id = abs(hash("cross-section"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    # Adding plots
    title = "Density profiles (Central haloes)"
    caption = ""
    filename = "Density_profiles_halos.png"
    id = abs(hash("Density profiles central halos"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Density profiles (Subhaloes)"
    caption = ""
    filename = "Density_profiles_subhalos.png"
    id = abs(hash("Density profiles subhalos"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    # Setting section
    title = "Halo/Subhalo density profiles"
    id = abs(hash("Profiles-density"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

    title = "Cross-section profiles (Central haloes)"
    caption = ""
    filename = "Cross_section_profiles_halos.png"
    id = abs(hash("Cross section profiles central halos"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    title = "Cross-section profiles (Subhaloes)"
    caption = ""
    filename = "Cross_section_profiles_subhalos.png"
    id = abs(hash("Cross section profiles subhalos"))
    PlotsInWeb.load_plots(title, caption, filename, id)

    # Setting section
    title = "Halo/Subhalo cross-section profiles"
    id = abs(hash("Profiles-sigma"))
    plots = PlotsInWeb.plots_details
    caption = " "
    add_web_section(web, title, caption, id, plots)
    PlotsInWeb.reset_plots_list()

