function area_TLRR=AUC(E0,mask)
%     Show_Our = sqrt(sum(E0.^2,3));
%     Show_Our = (Show_Our-min(Show_Our(:)))/(max(Show_Our(:))-min(Show_Our(:)));
    Show_Our = E0;
    r_max = max(Show_Our(:)); num = size(mask,1)*size(mask,2);
    taus = linspace(0, r_max, 5000);
    mask_reshape = reshape(mask, 1, num);
    anomaly_map = logical(double(mask_reshape)>0);                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           
    normal_map = logical(double(mask_reshape)==0);
    for index2 = 1:length(taus)
        tau = taus(index2);
        anomaly_map_rx = (Show_Our(:)> tau)';
        PF0(index2) = sum(anomaly_map_rx & normal_map)/sum(normal_map);
        PD0(index2) = sum(anomaly_map_rx & anomaly_map)/sum(anomaly_map);
    end
    area_TLRR = sum((PF0(1:end-1)-PF0(2:end)).*(PD0(2:end)+PD0(1:end-1))/2);


end