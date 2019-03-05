function [qvzelj status] = primruke(qvzelj,flier,brz,i) %#ok<INUSL>

    global rkukl rkukd rpetal rpetad rramel rramed rrukal rrukad 
    global qnew seglruka segdruka r_saka_d_p r_saka_l_p 
    
% %     r_saka_d_p = rrukad;
% %     r_saka_l_p = rrukal;
% % 
% %     rz = 0.69;
% % 
% %     r_saka_d_p(1) = rramed(1) + rskzl(1) - rkukl(1)+0.0783;% r_saka_rame_d_e;
% %     r_saka_l_p(1) = rramel(1) + rskzd(1) - rkukd(1)+0.0783;% r_saka_rame_l + r_saka_rame_l_e;

    r_saka_d_p = rramed;
    r_saka_l_p = rramel;
    
    rpetal_kor = rpetal(1:3) + RotX(rpetal(4))*RotY(rpetal(5))*RotZ(rpetal(6))*[0.1686;-0.28;0];
    rpetad_kor = rpetad(1:3) + RotX(rpetad(4))*RotY(rpetad(5))*RotZ(rpetad(6))*[0.1686;0.28;0];
    
    r_saka_d_p(1) = rpetal_kor(1);
    r_saka_l_p(1) = rpetad_kor(1);

    r_saka_d_p(2) = rramed(2);
    r_saka_l_p(2) = rramel(2);

    red = r_saka_d_p - rramed;
    rel = r_saka_l_p - rramel;
    
    red2 = red.*red;
    rel2 = rel.*rel;
    
    r_saka_d_p(3) = rramed(3)-sqrt(0.76*0.76 -red2(1)-red2(2));
    r_saka_l_p(3) = rramel(3)-sqrt(0.76*0.76 -rel2(1)-rel2(2));

    
    
%     r_saka_d_p(1) = rramed(1) + 1*(rpetal_kor(1) - rkukl(1));% r_saka_rame_d_e;
%     r_saka_l_p(1) = rramel(1) + 1*(rpetad_kor(1) - rkukd(1));% r_saka_rame_l + r_saka_rame_l_e;
% 
%     
    r_ruka_d_e = r_saka_d_p - rrukad;
    r_ruka_l_e = r_saka_l_p - rrukal;
    
    r_ruka_d_ort = [(r_ruka_d_e(1:3))/VecMod(r_ruka_d_e(1:3)) ; (r_ruka_d_e(4:6))/VecMod(r_ruka_d_e(4:6))];                                       
    r_ruka_l_ort = [(r_ruka_l_e(1:3))/VecMod(r_ruka_l_e(1:3)) ; (r_ruka_l_e(4:6))/VecMod(r_ruka_l_e(4:6))];                                       
    
    r_ruka_d_ort(isnan(r_ruka_d_ort))=0;
    r_ruka_l_ort(isnan(r_ruka_l_ort))=0;    
        
    if(brz>300*VecMod(r_ruka_d_e(1:3)))
        brzrd=300*VecMod(r_ruka_d_e(1:3));
    else
        brzrd=brz;
    end
        
    if(brz>300*VecMod(r_ruka_l_e(1:3)))
        brzrl=300*VecMod(r_ruka_l_e(1:3));
    else
        brzrl=brz;
    end
    
    qvrd = brzrd * r_ruka_d_ort;
    qvrl = brzrl * r_ruka_l_ort;
    
    [Jkrukad Akrukad] =km_jakP(flier,22);
    qvzelj(segdruka([1 3:4]),i)= pinv(Jkrukad(1:3,segdruka([1 3:4])))*qvrd(1:3);

    if(qnew(segdruka(4),i)>deg2rad(-10))
        qvzelj(segdruka(4),i)=0;
    end
    
    [Jkrukal Akrukal] =km_jakP(flier,23);
    qvzelj(seglruka([1 3:4]),i)= pinv(Jkrukal(1:3,seglruka([1 3:4])))*qvrl(1:3);
    
    if(qnew(seglruka(4),i)>deg2rad(-10))
        qvzelj(seglruka(4),i)=0;
    end
    

    
    
    
% 
%     
%     
%     K1 = -1;
%     K2 = -2;
%     K3 = -2;
%     
%     if(r_saka_d_p(1)>0)
%         q_z_rd = K1*r_saka_d_p(1);
%         q_z_ld = 0;
%     else
%         q_z_rd = K2*r_saka_d_p(1);
%         q_z_ld = K3*r_saka_d_p(1);
%     end
%         
%     if(r_saka_l_p(1)>0)
%         q_z_rl = K1*r_saka_l_p(1);
%         q_z_ll = 0;K2*r_saka_l_p(1);
%     else
%         q_z_rl = K2*r_saka_l_p(1);
%         q_z_ll = K3*r_saka_d_p(1);
%     end
%     
%     e_q_z_rd = q_z_rd - qnew(segdruka(3),i);
%     e_q_z_ld = q_z_ld - qnew(segdruka(4),i);
%     
%     e_q_z_rl = q_z_rl - qnew(seglruka(3),i);
%     e_q_z_ll = q_z_ll - qnew(seglruka(4),i);
%     
%     
%     Kp=-1;
%     
%     qvzelj(segdruka(3),i)=Kp*e_q_z_rd;
%     qvzelj(segdruka(4),i)=Kp*e_q_z_ld;
%     
%     qvzelj(seglruka(3),i)=Kp*e_q_z_rl;
%     qvzelj(seglruka(4),i)=Kp*e_q_z_ll;
    status=1;
    
end