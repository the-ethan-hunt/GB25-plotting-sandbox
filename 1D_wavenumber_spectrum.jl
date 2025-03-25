using Polynomials

function detrend_linear(x)
	p = fit(Polynomial, 1:length(x), x, 1)  
        trend = p(1:length(x))  
        return x .- trend
     end

using NCDatasets, FFTW, Statistics, Plots, DSP
     
function function_1d_spectrum(filename, latitude)
                  latitude_bounds = (-75,-25)
                  longitude_bounds = (200,250)
                  Nx = 50*4
                  longitude_span = longitude_bounds[2] - longitude_bounds[1]
                  dxF_deg = longitude_span / Nx  # in degrees
                  dxF_km = dxF_deg * cos(latitude_bounds[1] * π / 180) * 111
                  ds = Dataset(filename)
                  v = ds["v"][:, :,:, end]
                  yF = ds["yF"][:]
                  yF_index = argmin(abs.(yF .- latitude))
                  v_selected = v[:, yF_index, : ,:]
                  v_selected = replace(v_selected, NaN => 0.0)
                  v_mean_time = mean(v, dims=1)
                  v_centered = v .- v_mean_time
                  v_detrended = detrend_linear(v_centered)
                  v_ft = fft(v_detrended, dims=2)
                  viso2 = abs.(v_ft).^2
                  window = hanning(size(v, 2))
                  v_windowed = v .* window'
                  v_ft_windowed = fft(v_windowed, dims=2)
                  viso2_windowed = abs.(v_ft_windowed).^2
                  ekeiso = 0.5 * viso2
                  nk = div(length(ekeiso[1, :]), 2)
                  ekeiso_half = ekeiso[:, 1:nk]
                  wavenumber = (1:nk) * 1000 * 2 * π
                  ekeiso_mean = mean(ekeiso, dims=1)
                  ekeiso_std = std(ekeiso, dims=1)
                  ekeiso_mean = squeeze(ekeiso_mean)
                  ekeiso_std = squeeze(ekeiso_std)
                  plot(wavenumber, ekeiso_mean, linewidth=3, label="Time mean", color=:black)
                  plot!(wavenumber, ekeiso_mean .- ekeiso_std, fillrange=ekeiso_mean .+ ekeiso_std, fillalpha=0.3, label="Envelope (±σ)", color=:orange)
                  ref_k = 20
                  k_ref = wavenumber[ref_k]
                  E_ref = ekeiso_mean[ref_k]
                  k_slope = LinRange(wavenumber[5], wavenumber[end], 10)
                  E_slope = E_ref .* (k_slope ./ k_ref) .^ -3
                  plot(k_slope, E_slope, linestyle=:dash, color=:green, linewidth=2, label="k^{-3} slope")
                  savefig("k_slope_fig.png")
                  end

filename = "target_dir/quarter_degree_velocity_averages.nc"

latitude = -60                  

function_1d_spectrum(filename, latitude)
