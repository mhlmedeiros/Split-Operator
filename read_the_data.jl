using HDF5;




primeira_camada = h5open("dados_split_operator.h5","r") do file
    names(file)
end






segunda_camada = h5open("dados_split_operator.h5","r") do file
    segunda_camada = [];
    for name in primeira_camada
        append!(segunda_camada,names(file[name]))
    end
    return segunda_camada
end





terceira_camada = h5open("dados_split_operator.h5","r") do file
    terceira_camada = [];
    for name in primeira_camada
        g = file[name];
        for name_grupo in segunda_camada
            append!(terceira_camada,names(g[name_grupo]))
        end
    end
    return terceira_camada
end





V_lido = h5open("dados_split_operator.h5","r") do file
    group = file[primeira_camada[1]];
    subgroup = group[segunda_camada[1]];
    read(subgroup["Potencial"])
end





Atr_grupo = h5readattr("dados_split_operator.h5", "split_1D/potencial_degrau")
#
#
#
# h5writeattr("dados_split_operator.h5", "split_1D/potencial_degrau",
#         Dict( "x_0" => string(round(x₀,1)),
#               "k_0" => string(round(k₀,1)),
#               "sigma" => string(round(σ,1))));
#
h5readattr("dados_split_operator.h5","split_1D/potencial_degrau/Potencial")
#
# name_path_real = string("split_1D/potencial_degrau/",nome_dados_real)
# name_path_imag = string("split_1D/potencial_degrau/",nome_dados_imag)
#
# psi_lido_real = h5read("dados_split_operator.h5",name_path_real)
# psi_lido_imag = h5read("dados_split_operator.h5",name_path_imag)
#
# psi_lido_total = complex(psi_lido_real,psi_lido_imag)



# plot(x_grid,abs(ψ_1).^2)
# plot(x_grid,0.5*10.0^(-2)*V_x, linewidth=2.0);
# # ylim(-0.05,1.05);
# # xlim(-20,20)
# title("Gaussian distribution");
# xlabel(L"$x$");
# ylabel(L"$|\psi(x)|^2$");
