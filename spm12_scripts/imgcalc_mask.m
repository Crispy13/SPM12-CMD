%-----------------------------------------------------------------------
% Job saved on 20-May-2020 14:07:53 by cfg_util (rev $Rev: 7345 $)
% spm SPM - Unknown
% cfg_basicio BasicIO - Unknown
%-----------------------------------------------------------------------
matlabbatch{1}.spm.util.imcalc.input = {
                                        '/data/eck/T1CE-T2_images/15782676/2016-07-12/T1CE_AX_1mm/805_t1axgd_1m.nii,1'
                                        '/data/eck/T1CE-T2_images/15782676/2016-07-12/T1CE_AX_1mm/mean805_t1axgd_1m.nii,1'
                                        '/data/eck/T1CE-T2_images/15782676/2016-07-12/T1CE_AX_1mm/r805 T1AXGD 1M_roi_label.nii,1'
                                        '/data/eck/T1CE-T2_images/15782676/2016-07-12/T1CE_AX_1mm/805 T1AXGD 1M_roi_label.nii,1'
                                        };
matlabbatch{1}.spm.util.imcalc.output = 'output';
matlabbatch{1}.spm.util.imcalc.outdir = {'/data/eck/T1CE-T2_images/15782676/2016-07-12/T1CE_AX_1mm'};
matlabbatch{1}.spm.util.imcalc.expression = '((i2+i3+i4)>0.5)';
matlabbatch{1}.spm.util.imcalc.var = struct('name', {}, 'value', {});
matlabbatch{1}.spm.util.imcalc.options.dmtx = 0;
matlabbatch{1}.spm.util.imcalc.options.mask = 0;
matlabbatch{1}.spm.util.imcalc.options.interp = 1;
matlabbatch{1}.spm.util.imcalc.options.dtype = 4;
