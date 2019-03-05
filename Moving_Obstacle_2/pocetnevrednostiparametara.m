function pocetnevrednostiparametara

    %Prvi polukorak
    %faza 1
    Kb_osl_L = 0.3; %#ok<NASGU>
    Kv_osl_L = 0.93; %#ok<NASGU>

    %faza 1.5
    Kb_nn_L = 0.075; %#ok<NASGU>
    Kv_nn_L = 0.93; %#ok<NASGU>
    
    ikrajnn_L = 0; %#ok<NASGU>
    
    Kb_gn_D = 0.6; %#ok<NASGU>
    Kv_gn_D = 0.1; %#ok<NASGU>
    Kd_gn_D = 0.02; %#ok<NASGU>

    ikrajgd_D = 0; %#ok<NASGU>

    %faza 2    
    Kb_nn2_L = 0.1; %#ok<NASGU>
    Kv_nn2_L = 0.93; %#ok<NASGU>
    
    ikrajnn2_L = 0; %#ok<NASGU>
    
    Kb_on_D = 0.3; %#ok<NASGU>
    Kv_on_D = 0.05; %#ok<NASGU>
    Kd_on_D = 0.20; %#ok<NASGU>
    Ku_on_D = -0.1; %#ok<NASGU>
        
    ikrajon_D = 0; %#ok<NASGU>

    %faza 2.5
    Kb_nn25_L = 0.1; %#ok<NASGU>
    Kv_nn25_L = 0.93; %#ok<NASGU>
    
    ikrajnn25_L = 0; %#ok<NASGU> 
    
    Kb_on25_D = 0.4; %#ok<NASGU>
    Kv_on25_D = 0.005; %#ok<NASGU>
    Kd_on25_D = 0.22; %#ok<NASGU>
    Ku_on25_D = -0.08; %#ok<NASGU>

    ikrajon25_D = 0; %#ok<NASGU>
    
    %faza 3
    Kb_sn_D = 0.3; %#ok<NASGU>
    Kv_sn_D = 0.93; %#ok<NASGU>
    
     %Drugo polukorak
    %faza 4
    Kb_osl_D = 0.3; %#ok<NASGU>
    Kv_osl_D = 0.93; %#ok<NASGU>

    %faza 4.5
    Kb_nn_D = 0.075; %#ok<NASGU>
    Kv_nn_D = 0.93; %#ok<NASGU>
    
    ikrajnn_D = 0; %#ok<NASGU>
    
    Kb_gn_L = 0.6; %#ok<NASGU>
    Kv_gn_L = 0.1; %#ok<NASGU>
    Kd_gn_L = 0.02; %#ok<NASGU>

    ikrajgn_L = 0; %#ok<NASGU>

    %faza 5    
    Kb_nn5_D = 0.1; %#ok<NASGU>
    Kv_nn5_D = 0.93; %#ok<NASGU>
    
    ikrajnn5_D = 0; %#ok<NASGU>

    Kb_on_L = 0.3; %#ok<NASGU>
    Kv_on_L = 0.05; %#ok<NASGU>
    Kd_on_L = 0.20; %#ok<NASGU>
    Ku_on_L = -0.1; %#ok<NASGU>
        
    ikrajon_L = 0; %#ok<NASGU>

    %faza 5.5
    Kb_nn55_D = 0.1; %#ok<NASGU>
    Kv_nn55_D = 0.93; %#ok<NASGU>
    
    ikrajnn55_D = 0; %#ok<NASGU>

    Kb_on55_L = 0.4; %#ok<NASGU>
    Kv_on55_L = 0.005; %#ok<NASGU>
    Kd_on55_L = 0.22; %#ok<NASGU>
    Ku_on55_L = -0.08;  %#ok<NASGU>   
            
    ikrajon55_L = 0; %#ok<NASGU>

    %faza 6
    Kb_sn_L = 0.3; %#ok<NASGU>
    Kv_sn_L = 0.93; %#ok<NASGU>
    
    
    save pocetniparametri.mat
end