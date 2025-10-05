# Minimal, robust processing of the uploaded NASA Giovanni NetCDF (.nc).
# Handles both time series (with 'time' dimension) and time-averaged map (no 'time').
# Produces downloads when time series is present.

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Try importing xarray and netCDF4
xr = None
Dataset = None
try:
    import xarray as _xr
    xr = _xr
except Exception:
    xr = None

try:
    from netCDF4 import Dataset as _Dataset, num2date
    Dataset = _Dataset
except Exception:
    Dataset = None
    num2date = None

# (Apenas neste ambiente) para mostrar tabelas ao usuário.
# Se for rodar localmente, COMENTE a linha abaixo e as duas chamadas mais adiante.
from caas_jupyter_tools import display_dataframe_to_user

# >>> Troque o caminho abaixo para o seu arquivo .nc quando rodar localmente <<<
nc_path = "/mnt/data/g4.timeAvgMap.GPM_3IMERGM_07_precipitation.20250201-20251031.38W_11S_35W_8S.nc"

def open_dataset(nc_path):
    if xr is not None:
        try:
            return xr.open_dataset(nc_path), "xarray"
        except Exception:
            pass
    if Dataset is not None:
        ds = Dataset(nc_path, mode="r")
        return ds, "netcdf4"
    raise RuntimeError("Não consegui abrir o NetCDF com xarray ou netCDF4.")

def find_precip_var_xr(ds):
    for v in ds.data_vars:
        name = v
        var = ds[v]
        units = str(var.attrs.get("units","")).lower()
        long_name = str(var.attrs.get("long_name","")).lower()
        if "precip" in name.lower() or "precip" in long_name:
            return v
    # fallback
    return list(ds.data_vars.keys())[0]

def find_precip_var_nc(nc):
    # pick first variable that is not a coordinate and has lat/lon dims if possible
    cand = None
    for k, v in nc.variables.items():
        dims = getattr(v, "dimensions", ())
        if "lat" in [d.lower() for d in dims] or "latitude" in [d.lower() for d in dims]:
            cand = k
            if "precip" in k.lower():
                return k
    # weaker fallback: first non-dimension var
    for k, v in nc.variables.items():
        if k not in nc.dimensions:
            return k
    return None

def get_latlon_names(dims):
    lat_name = None
    lon_name = None
    for d in dims:
        dl = d.lower()
        if dl.startswith("lat"):
            lat_name = d
        if dl.startswith("lon"):
            lon_name = d
    return lat_name, lon_name

def coslat_weights(lat):
    w = np.cos(np.deg2rad(lat))
    w = np.clip(w, 0, None)
    s = np.nansum(w)
    return w / s if s != 0 else w

ds, backend = open_dataset(nc_path)

if backend == "xarray":
    dims = list(ds.dims)
    has_time = "time" in ds.dims
    lat_name, lon_name = get_latlon_names(ds.dims)
    if len(ds.data_vars) == 0:
        raise RuntimeError("Arquivo .nc não possui variáveis de dados.")
    var_name = find_precip_var_xr(ds)
    var = ds[var_name]
    units = var.attrs.get("units", "unknown")
    
    if has_time:
        # area-mean if lat/lon exist
        if lat_name in var.dims and lon_name in var.dims:
            lat = ds[lat_name].values
            w = coslat_weights(lat)
            W = np.repeat(w[:, None], ds[lon_name].size, axis=1)
            # weighted mean
            area_mean = var.weighted(xr.DataArray(W, dims=(lat_name, lon_name))).mean(dim=(lat_name, lon_name))
        else:
            area_mean = var

        # Convert to mm/month from mm/day
        s = pd.Series(area_mean.values, index=pd.to_datetime(area_mean["time"].values)).sort_index().dropna()
        s_mm_month = s * s.index.days_in_month

        # Climatology
        grp = s_mm_month.groupby(s_mm_month.index.month)
        clim = grp.mean()
        p10 = grp.quantile(0.10)
        p90 = grp.quantile(0.90)
        clim_df = pd.DataFrame({
            "climatology_mm_month": [clim.get(m, np.nan) for m in range(1,13)],
            "p10_mm_month": [p10.get(m, np.nan) for m in range(1,13)],
            "p90_mm_month": [p90.get(m, np.nan) for m in range(1,13)]
        }, index=pd.Index(range(1,13), name="month"))

        # Outputs
        out_csv = "/mnt/data/IMERG_area_mean_monthly.csv"
        out_xlsx = "/mnt/data/Historico_Climatico_from_NC.xlsx"
        out_png = "/mnt/data/Historico_Climatico_timeseries.png"

        s_mm_month.to_csv(out_csv, header=["mm_month"])

        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
            s_mm_month.to_frame("mm_month").to_excel(w, sheet_name="SERIE", index_label="date")
            clim_df.to_excel(w, sheet_name="HISTORICO_CLIMATOLOGIA")

        # Chart
        plt.figure(figsize=(10,5))
        plt.plot(s_mm_month.index, s_mm_month.values)
        plt.title("Histórico Climático — IMERG V07 (mm/mês)")
        plt.xlabel("Data")
        plt.ylabel("mm/mês")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out_png, dpi=140)
        plt.close()

        # Mostrar tabelas (REMOVA se rodar fora deste ambiente)
        display_dataframe_to_user("Série mensal (últimos 24 meses)", s_mm_month.tail(24).to_frame("mm_month").reset_index())
        display_dataframe_to_user("Climatologia por mês (média, P10, P90)", clim_df.reset_index())

        print("OK: Série temporal encontrada e processada.")
        print(f"[Download CSV] {out_csv}")
        print(f"[Download Excel] {out_xlsx}")
        print(f"[Download PNG] {out_png}")

    else:
        # time-averaged map (no time) — cannot make line chart
        lat_name, lon_name = get_latlon_names(var.dims)
        if lat_name and lon_name:
            mean_val = float(var.mean().values)
        else:
            mean_val = float(np.nanmean(var.values))
        out_csv = "/mnt/data/IMERG_timeAvgMap_mean.csv"
        pd.DataFrame([{"metric":"area_mean_precip", "value":mean_val, "units":units}]).to_csv(out_csv, index=False)
        print("ATENÇÃO: O .nc é um timeAvgMap (média no tempo); não possui dimensão temporal para série histórica.")
        print(f"Média regional do mapa: {mean_val:.3f} {units}")
        print(f"[Download CSV] {out_csv}")

else:
    # netCDF4 backend processing
    vars_list = list(ds.variables.keys())
    dims_list = list(ds.dimensions.keys())
    has_time = "time" in ds.dimensions
    var_name = find_precip_var_nc(ds)
    var = ds.variables[var_name]
    units = getattr(var, "units", "unknown")
    # Find lat/lon dim names
    lat_name, lon_name = None, None
    for d in var.dimensions:
        dl = d.lower()
        if dl.startswith("lat"):
            lat_name = d
        if dl.startswith("lon"):
            lon_name = d

    if has_time:
        arr = var[:]
        if lat_name and lon_name and lat_name in var.dimensions and lon_name in var.dimensions:
            lat = ds.variables[lat_name][:]
            w = coslat_weights(lat)
            W = np.repeat(w[:, None], ds.dimensions[lon_name].size, axis=1)
            # Weighted mean across lat/lon
            dims = var.dimensions
            t_idx = dims.index("time")
            lat_idx = dims.index(lat_name) if lat_name in dims else None
            lon_idx = dims.index(lon_name) if lon_name in dims else None

            # Reorder to (time, lat, lon) if needed
            perm = [t_idx, lat_idx, lon_idx]
            if len(dims) == 3 and all(i is not None for i in perm):
                arr_tll = np.moveaxis(arr, perm, [0,1,2])
                weighted = arr_tll * W[None, :, :]
                area_mean_vals = np.nanmean(weighted, axis=(1,2))
            else:
                area_mean_vals = np.nanmean(arr, axis=tuple(i for i in range(arr.ndim) if i!=t_idx))
        else:
            # Already a time series
            dims = var.dimensions
            t_idx = dims.index("time")
            area_mean_vals = np.nanmean(var[:], axis=tuple(i for i in range(var[:].ndim) if i!=t_idx))

        # build time index
        time_var = ds.variables["time"]
        if num2date is not None and hasattr(time_var, "units"):
            time_vals = num2date(time_var[:], units=time_var.units)
            time_index = pd.to_datetime([pd.Timestamp(t).to_pydatetime() for t in time_vals])
        else:
            time_index = pd.to_datetime(time_var[:], errors="coerce")

        s_mm_day = pd.Series(area_mean_vals, index=time_index).sort_index().dropna()
        s_mm_month = s_mm_day * s_mm_day.index.days_in_month

        grp = s_mm_month.groupby(s_mm_month.index.month)
        clim = grp.mean()
        p10 = grp.quantile(0.10)
        p90 = grp.quantile(0.90)
        clim_df = pd.DataFrame({
            "climatology_mm_month": [clim.get(m, np.nan) for m in range(1,13)],
            "p10_mm_month": [p10.get(m, np.nan) for m in range(1,13)],
            "p90_mm_month": [p90.get(m, np.nan) for m in range(1,13)]
        }, index=pd.Index(range(1,13), name="month"))

        out_csv = "/mnt/data/IMERG_area_mean_monthly.csv"
        out_xlsx = "/mnt/data/Historico_Climatico_from_NC.xlsx"
        out_png = "/mnt/data/Historico_Climatico_timeseries.png"

        s_mm_month.to_csv(out_csv, header=["mm_month"])

        with pd.ExcelWriter(out_xlsx, engine="openpyxl") as w:
            s_mm_month.to_frame("mm_month").to_excel(w, sheet_name="SERIE", index_label="date")
            clim_df.to_excel(w, sheet_name="HISTORICO_CLIMATOLOGIA")

        plt.figure(figsize=(10,5))
        plt.plot(s_mm_month.index, s_mm_month.values)
        plt.title("Histórico Climático — IMERG V07 (mm/mês)")
        plt.xlabel("Data")
        plt.ylabel("mm/mês")
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(out_png, dpi=140)
        plt.close()

        # Mostrar tabelas (REMOVA se rodar fora deste ambiente)
        display_dataframe_to_user("Série mensal (últimos 24 meses)", s_mm_month.tail(24).to_frame("mm_month").reset_index())
        display_dataframe_to_user("Climatologia por mês (média, P10, P90)", clim_df.reset_index())

        print("OK: Série temporal encontrada e processada (backend netCDF4).")
        print(f"[Download CSV] {out_csv}")
        print(f"[Download Excel] {out_xlsx}")
        print(f"[Download PNG] {out_png}")

    else:
        # time averaged map — no time dim
        data = var[:]
        mean_val = float(np.nanmean(data))
        out_csv = "/mnt/data/IMERG_timeAvgMap_mean.csv"
        pd.DataFrame([{"metric":"area_mean_precip", "value":mean_val, "units":units}]).to_csv(out_csv, index=False)
        print("ATENÇÃO: O .nc é um timeAvgMap (média no tempo); não possui dimensão temporal para série histórica.")
        print(f"Média regional do mapa: {mean_val:.3f} {units}")
        print(f"[Download CSV] {out_csv}")

