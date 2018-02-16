using HDF5;

function save_all(h, g, subg_name, ψ, data_name, ini_cond_attr,  pot_attr)
    save_potential(h, pot_attr)
    save_cond_ini(h, g, subg_name, ini_cond_attr)
    save_function(h, ψ, data_name)
end

function save_potential(h, attr_dict)
    h["Potencial"] = V_x;
    for (key, value) in attr_dict
        # Os atributos do Potencial são dados pelas tuplas (key, value) do
        # dictionary "attr_dict".
        attrs(h["Potencial"])[key] = string(value);
    end
end

function save_cond_ini(h, g, subg_name, attr_dict)
    # Salva a condição inicial:
    h["Cond_Ini_real"] = real(Ψ_0);
    h["Cond_Ini_imag"] = imag(Ψ_0);
    h["grid"] = collect(x_grid);
    for (key, value) in attr_dict
        attrs(g[subg_name])[key] = string(value);
    end
end

function save_function(handle, ψ, data_name)
    # Salva a solução para o tempo final:
    nome_real = string(data_name,"_real");
    nome_imag = string(data_name,"_imag");
    handle[nome_real] = real(ψ);
    handle[nome_imag] = imag(ψ);
end

function save_them_all(file_name, group_name, subg_name, ψ, data_name, ini_cond_attr,  pot_attr)
    try
        h5open(file_name,"r+") do file
            # Cria os grupos, ou abre aqueles já existentes:
            try
                g = file[group_name];
                h = g[ subg_name];
                save_function(h, ψ, data_name)
            catch
                g = g_create(file, group_name);
                h = g_create(g, subg_name);
                save_all(h, g, subg_name, ψ, data_name, ini_cond_attr,  pot_attr)
            end
        end
    catch
        h5open(file_name,"w") do file
            # Cria os grupos
            g = g_create(file, group_name);
            h = g_create(g, subg_name);
            save_all(h, g, subg_name, ψ, data_name, ini_cond_attr,  pot_attr)
        end
    end
end
