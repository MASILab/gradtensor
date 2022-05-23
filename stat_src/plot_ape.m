ape = nifti_utils.load_untouch_nii4D_vol_scaled('/nfs/masi/kanakap/projects/LR/masivar/fa_ape_image.nii','double');
slice = ape(48,:,:);
figure; imshow(rot90(squeeze(slice)),[0, 10], 'Colormap', hot); colorbar;
ax = gca;
ax.FontSize = 17;
