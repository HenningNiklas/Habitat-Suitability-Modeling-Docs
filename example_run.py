from habitatsuitability import run_complete_hsi_habitat_binning
import os



def main():
    velocity_hsicurve_json_path = os.path.abspath("") + "\\example_inputs\\exampleFish_hsi_curve_json\\_example_velocity_hsi.json"
    depth_hsicurve_json_path = os.path.abspath("") + "\\example_inputs\\exampleFish_hsi_curve_json\\_example_depth_hsi.json"
    substrate_hsicurve_json_path = os.path.abspath("") + "\\example_inputs\\exampleFish_hsi_curve_json\\_example_substrate_hsi.json"
    other_hsicurve_json_path = None
    parameter_hsicurve_json_dict = dict(
        velocity=velocity_hsicurve_json_path,
        depth=depth_hsicurve_json_path,
        substrate=substrate_hsicurve_json_path,
        other=other_hsicurve_json_path)
    tif_inputs = {"velocity": os.path.abspath("") + "\\example_inputs\\input_parameter_tifs\\example_velocity.tif",
                  "depth": os.path.abspath("") + "\\example_inputs\\input_parameter_tifs\\example_depth.tif",
                  "substrate": os.path.abspath("") + "\\example_inputs\\input_parameter_tifs\\example_substrate.tif"}

    run_complete_hsi_habitat_binning(tifs_dictionary=tif_inputs,  # Input test +create_hsi_raster
                                     parameters=["velocity", "depth", "substrate"],  # create_hsi_raster
                                     parameter_hsicurve_json_dict=parameter_hsicurve_json_dict,
                                     # bioverification
                                     depth_tif_path=tif_inputs["depth"],
                                     observation_coordinate_csv_path=os.path.abspath(
                                         "") + "\\example_inputs\\input_bioverification_files\\example_observation_data.csv",

                                     bin_resolution=0.04, random_datasets_quantity=10000)  # habitat_binning)




    #print(name)
if __name__ == '__main__':


    main()

