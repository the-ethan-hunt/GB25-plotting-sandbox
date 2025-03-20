import numpy as np
import xarray as xr

def latitude_longitude_grid(ds, resolution):

    # Earth's radius in meters (approximate mean radius)
    R = 6371.0 * 1000  # Convert kilometers to meters
    
    # Assume ds["yC"] and ds["xC"] represent latitude and longitude arrays in degrees
    lat = ds["yC"]
    latF = ds["yF"]
    lon = ds["xC"]
    
    # Create an empty dataset
    ds_grid = xr.Dataset()
    
    # Convert degrees to radians
    lat_rad = np.deg2rad(lat)
    latF_rad = np.deg2rad(latF)
    lon_rad = np.deg2rad(lon)
    
    # Compute y-spacing (difference in latitude)
    # Each degree of latitude corresponds to approximately the same distance
    y_spacing = (np.deg2rad(resolution) * R).item()  # Convert degree of latitude to meters
    
    # Compute x-spacing (difference in longitude)
    # Distance depends on the latitude
    x_spacing = R * np.deg2rad(resolution) * np.cos(lat_rad)
    xF_spacing = R * np.deg2rad(resolution) * np.cos(latF_rad)
    # Add spacing to the dataset
    ds_grid["dy"] = xr.DataArray(y_spacing, dims=[])
    ds_grid["dy"].attrs["long_name"] = "Grid spacing in y-direction"
    ds_grid["dy"].attrs["units"] = "m"

    ds_grid["dx"] = xr.DataArray(x_spacing, dims=["yC"], coords={"yC": lat})
    ds_grid["dx"].attrs["long_name"] = "Grid spacing in x-direction"
    ds_grid["dx"].attrs["units"] = "m"
    
    ds_grid["dxF"] = xr.DataArray(xF_spacing, dims=["yF"], coords={"yF": latF})
    ds_grid["dxF"].attrs["long_name"] = "Grid spacing in x-direction"
    ds_grid["dxF"].attrs["units"] = "m"
    
    dz = xr.DataArray(np.diff(ds["zF"]), dims=["zC"])
    ds_grid["dz"] = dz
    ds_grid["dz"].attrs["long_name"] = "Grid spacing in z-direction"
    ds_grid["dz"].attrs["units"] = "m"
    
    if "eta" in ds.data_vars:
        eta = xr.zeros_like(dz) * xr.zeros_like(ds.eta)
        eta.loc[{"zC": eta.zC[-1].item()}] = ds.eta
        ds_grid["dz_adjusted"] = dz + eta

    # Define a helper function for `hack_sind`
    def hack_sind(degrees):
        return np.sin(np.deg2rad(degrees))
    
    # Compute area in square meters
    ds_grid["area"] = (
        R**2 *
        np.deg2rad(resolution) * hack_sind(ds["yF"]).diff(dim="yF").swap_dims({"yF": "yC"})
    ).assign_coords({"yC": ds.yC})
    ds_grid["area"].attrs["long_name"] = "Grid cell area"
    ds_grid["area"].attrs["units"] = "m2"

    return ds_grid

