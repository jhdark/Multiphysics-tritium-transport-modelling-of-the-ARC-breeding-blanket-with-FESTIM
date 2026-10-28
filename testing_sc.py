from festim_script import build_festim_model, load_openfoam_data, define_meshtags

if __name__ == "__main__":
    run_filename = "data/NewOF_symmBCs/arc.foam"
    run_time = 745

    # load data
    mesh, u, T, nut, t_prod, cfd_reader = load_openfoam_data(
        filename=run_filename, time=run_time
    )
    # convert tritprod from kg/m3/s to T/m3/s
    tritium_at_mass = 3.01604928  # Da or u
    tritium_at_mass = tritium_at_mass * 1.660539e-27  # kg
    t_prod.x.array[:] /= tritium_at_mass  # T/m3/s

    # define meshtags
    facet_meshtags, volume_meshtags = define_meshtags(cfd_reader)

    # for Sc_value in [0.3]:
    for Sc_value in [0.3, 0.5, 0.7, 0.9, 1.1, 1.3]:
        print(f"Running with Turbuelnt schmidt number = {Sc_value}")
        my_model = build_festim_model(
            mesh=mesh,
            facet_meshtags=facet_meshtags,
            volume_meshtags=volume_meshtags,
            u=u,
            T=T,
            nut=nut,
            tritium_production_rate=t_prod,
            delta=0.1,
            Sc_no=Sc_value,
            export_foldername=f"results/symmBCs/sc_testing/sc_{Sc_value}",
            steady=True,
        )

        print("Initialising and running festim model")
        my_model.initialise()
        my_model.run()
