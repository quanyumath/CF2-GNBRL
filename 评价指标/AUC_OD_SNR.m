function [AUC_OD_TLRR, AUC_SNR_TLRR] = AUC_OD_SNR(E, mask)

if size(E, 3) > 1
    E = sqrt(sum(E.^2, 3));
end
E = (E - min(E(:))) / (max(E(:)) - min(E(:)));

[PF_TLRR, PD_TLRR, Tau_TLRR] = perfcurve(mask(:), E(:), '1');
AUC_TLRR = -sum((PF_TLRR(1:end-1) - PF_TLRR(2:end)).*(PD_TLRR(2:end) + PD_TLRR(1:end-1))/2);
AUC_TLRR_PDtau = sum((Tau_TLRR(1:end-1) - Tau_TLRR(2:end)).*(PD_TLRR(2:end) + PD_TLRR(1:end-1))/2);
AUC_TLRR_PFtau = sum((Tau_TLRR(1:end-1) - Tau_TLRR(2:end)).*(PF_TLRR(2:end) + PF_TLRR(1:end-1))/2);
AUC_OD_TLRR = AUC_TLRR + AUC_TLRR_PDtau - AUC_TLRR_PFtau;
AUC_SNR_TLRR = AUC_TLRR_PDtau / AUC_TLRR_PFtau;

end