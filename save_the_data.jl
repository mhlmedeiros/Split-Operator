using HDF5;

function save_data(name_file, name_group, name_sub_group)

    # Criar as coisas: Arquivo, Grupo, SubGrupo

    h5open(name_file,"w") do file
        # Cria, ou abre o grupo e o subgrupo:
        try g = g_create(file, name_group) catch g = g_open(file, name_group) end;
        try h = g_create(g, name_sub_group) catch h = g_open(file, name_sub_group) end;

        # Caso não tenha salvo o potencial, salve-o:
        try
            h["Potencial"] = V_x;
            attrs(h["Potencial"])["alpha"] = string(α);
            attrs(h["Potencial"])["beta"] = string(β);
            attrs(h["Potencial"])["gamma"] = string(γ);
        end
