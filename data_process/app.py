import xarray as xr
import pandas as pd

file_path = "df/SMAP_L3_SM_P_E_20251003_R19240_001.h5"
ds = xr.open_dataset(file_path, engine="h5netcdf", group="Soil_Moisture_Retrieval_Data_PM", decode_times=False, phony_dims="sort")

sm = ds["soil_moisture_pm"].values.astype(float)
lat = ds["latitude_pm"].values.astype(float)
lon = ds["longitude_pm"].values.astype(float)

# Visualização de variáveis' nomes
print(ds.variables)

# Gera CSV com os dados extraídos para análise exploratória, limpeza e retorno de insights
data = {
    "Longitude": lon.flatten(),
    "Latitude": lat.flatten(),
    "Soil_Moisture": sm.flatten()
}
df = pd.DataFrame(data)
df.to_csv('df/SMAP_L3_SM_P_E_20251003_R19240_001.csv', index=False)