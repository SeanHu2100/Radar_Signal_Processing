% plot range-angle heatmap
% function [axh] = plot_rangeDop(Dopdata_sum,rng_grid,vel_grid)
% 
% % plot 2D(range-Doppler)
% figure('visible','on')
% set(gcf,'Position',[10,10,530,420])
% [axh] = surf(vel_grid,rng_grid,Dopdata_sum);
% view(0,90)
% axis([-8 8 0.5 5.5]);
% grid off
% shading interp
% xlabel('Doppler Velocity (m/s)')
% ylabel('Range(meters)')
% colorbar
% % caxis([0,3e04])
% title('Range-Doppler heatmap')
% 
% end

function [axh] = plot_rangeDop(Dopdata_sum, rng_grid, vel_grid, Resl_indx)

% plot 2D(range-Doppler)
figure('visible','on');
set(gcf, 'Position', [10,10,530,420]);
axh = surf(vel_grid, rng_grid, Dopdata_sum);
view(0,90);
axis([-12 12 0 6]);
grid off;
shading interp;
xlabel('Doppler Velocity (m/s)');
ylabel('Range (meters)');
colorbar;
title('Range-Doppler Heatmap');

% ==== 可选：叠加 CFAR 检测目标点 ====
if nargin >= 4 && ~isempty(Resl_indx)
    hold on;
    for k = 1:size(Resl_indx, 2)
        dop_idx = Resl_indx(1, k);
        rng_idx = Resl_indx(2, k);
        
        % 映射到物理坐标
        vel = vel_grid(dop_idx);
        rng = rng_grid(rng_idx);
        
        % 红色圆点
        plot3(vel, rng, max(Dopdata_sum(:)), 'ro', 'MarkerSize', 8, 'LineWidth', 1.5);
    end
    hold off;
end

end
