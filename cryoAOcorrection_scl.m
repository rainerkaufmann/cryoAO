% cryoAOcorrection_scl(mmc,dm,Z2C,TH,Ws,Wa,pri_only,iter,summode,blcorr,savename)
%
%
% Copyright 2024-2026 by Rainer Kaufmann
% All rights reserved.
%
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%


function [Zsum,pbest,pall,img_start,img_end] = cryoAOcorrection_scl(mmc,dm,Z2C,TH,Wi,Ws,Wa,pri_only,prifirst,iter,summode,blcorr,savename)

% mirrorSN = 'BAX564';
% dm = asdkDM( mirrorSN );


% reduce Z2C to modes actually used
Z2C_max = Z2C(1:21,:);
if pri_only==1
    Zmax = 7;
else
    Zmax = 18;
end

mmc.setProperty("Andor", "Region of Interest","2. 256 x 256 (centered)")
% mmc.setProperty("Andor", "Region of Interest","3. 128 x 128 (centered)")
% mmc.setProperty("Andor", "Region of Interest","4. 64 x 64 (centered)")

dm.Reset();
pause(0.02)


% img = MM_getimg(mmc); % 1st image to clear memory
% pause(0.01) % wait for data capture
img = MM_getimg(mmc);
img_start = img;
CRstack = zeros([size(img),3],'uint16');
pause(0.02)
img_start = img;
figure('Position',[500 100 800 800])
imhandle = imagesc(img);
axis equal
colorbar

p = zeros(21,1)+0.0001;
pall = p;
p_i = p;
% a-index:
%     1    2    3    4    5    6    7    8    9    10   11   12   13   14
a = [0.00 0.00 0.00 0.30 0.30 0.50 0.50 0.30 0.30 0.50 0.15 0.15 0.15 0.15...
    0.25 0.25 0.15 0.15 0.15 0.15 0.15];
%    15   16   17   18   19   20   21

% Zernike mode order for correction routine
ZMO = [10 6 7 4 5 8 9 21 15 16 11 12 17 18 13 14 19 20];

Zyticklabels = {'1','2','(3) defocus','(4) asti_o_b_l','(5) asti_v_e_r',...
    '(6) coma_v_e_r','(7) coma_h_o_r','(8) tref_v_e_r','(9) tref_h_o_r',...
    '(10) spher', '(11) sec asti_v_e_r','(12) sec asti_o_b_l',...
    '(13) quad_v_e_r','(14) quad_o_b_l','(15) sec coma_h_o_r',...
    '(16) sec coma_v_e_r','(17) sec tref_o_b_l','(18) sec tref_v_e_r',...
    '(19) pent_o_b_l','(20) pent_v_e_r','(21) sec spher'};

if summode==1
    a = a.*0.6;
else
    a = a.*0.75;
end
start_a = a;

% asym limit
al = 2;


Ortes = POSfitter_vs(img_start,130,TH,5,1);
Ortel = POSfitter_vs(img_start,130,TH,9,1);
BeadsNum = size(Ortes,1)
Imid_best = mean(Ortel(:,11));
sig_best = mean(mean(Ortel(:,6:7)'));
asym_best = mean(Ortes(:,2:3) - Ortel(:,2:3),1);
% Qm_best = (Imid_best)./(Ws.*sig_best)./(Wa+sqrt((asym_best(1)).^2+(asym_best(2)).^2))
Qm_best = Wi*Imid_best/10000 + (Ws*1000)/sig_best + (Wa*10)/(al+sqrt((asym_best(1)).^2+(asym_best(2)).^2))
% Qm_best = (Imid_best)/sig_best;

Imid = Imid_best;
Imid_start = Imid;
sig = sig_best;
sig_start = sig;
asym = asym_best;
asym_start = mean(asym);
Qm = Qm_best;
Qm_start = Qm;
fighm = figure('Position',[120 120 800 1000]);
barhandle = subplot(4,2,[1,3,5,7]);
barh(barhandle,p(1:21));
xlim(barhandle,[-1 1])
ylim(barhandle,[3.5 21+0.5])
xlabel(barhandle,'amplitude (µm RMS)');
ylabel(barhandle,'Zernike mode');
yticks(barhandle,4:21)
yticklabels(barhandle,Zyticklabels(4:21));
figh1 = subplot(4,2,2);
iteration = 0;
plot(figh1,iteration,Imid,'o-');
yline(figh1,Imid_start,'k--');
xlabel(figh1,'iteration')
ylabel(figh1,'central intensity')
xlim(figh1,[0 iter])
figh2 = subplot(4,2,4);
plot(figh2,iteration,sig.*2.355,'o-');
yline(figh2,sig_start.*2.355,'k--');
xlabel(figh2,'iteration')
ylabel(figh2,'PSF width (nm)')
xlim(figh2,[0 iter])
figh3 = subplot(4,2,6);
plot(figh3,iteration,asym(2),'o-');
hold(figh3,'on')
plot(figh3,iteration,asym(1),'o-');
plot(figh3,iteration,sqrt(asym(1)^2+asym(2)^2),'o-');
hold(figh3,'off')
xlabel(figh3,'iteration')
yline(figh3,0,'k--');
ylabel(figh3,'PSF asymmetry (nm)')
% yline(figh3,asym_start,'k--');
legend(figh3,'x','y','xy','Location','NorthWest')
xlim(figh3,[0 iter])
figh4 = subplot(4,2,8);
plot(figh4,iteration,Qm,'o-');
yline(figh4,Qm_start,'k--');
xlabel(figh4,'iteration')
ylabel(figh4,'overall quality')
xlim(figh4,[0 iter])



for i=1:iter
    % aplly modes determined from previous iterations (0 at i=1) 
    p_s = p;
    Zsum = sum(Z2C_max.*p(1:21),1);
    for Zm=ZMO(1:Zmax)
        if Zm<11 | (Zm>10 & i>10) | prifirst==0
            % aplly for each mode amplitude -a, 0 and +a
    %         if summode==1
    %             dm.Send(Zsum + sum(Z2C_max.*p_i(1:21)));
    %         else
                dm.Send(DMscaling(Zsum));
    %         end
            pause(0.02)
            img = MM_getimg(mmc);
            CRstack(:,:,2) = img;
            set(imhandle, 'CData', img);
    %         if summode==1
    %             dm.Send(Zsum + Z2C(Zm,:).*a(Zm) + sum(Z2C_max.*p_i(1:21)));
    %         else
                dm.Send(DMscaling(Zsum + Z2C(Zm,:).*a(Zm)));
    %         end
            pause(0.02)
            img = MM_getimg(mmc);
            CRstack(:,:,3) = img;
            set(imhandle, 'CData', img);
    %         if summode==1
    %             dm.Send(Zsum - Z2C(Zm,:).*a(Zm) + sum(Z2C_max.*p_i(1:21)));
    %         else
                dm.Send(DMscaling(Zsum - Z2C(Zm,:).*a(Zm)));
    %         end
            pause(0.02)
            img = MM_getimg(mmc);
            CRstack(:,:,1) = img;
            set(imhandle, 'CData', img);

            % determine parameters of point-lilke-objects (PSFs)
            Ortes1 = POSfitter_vs(CRstack(:,:,1),130,TH,5,1);
            Ortes2 = POSfitter_vs(CRstack(:,:,2),130,TH,5,1);
            Ortes3 = POSfitter_vs(CRstack(:,:,3),130,TH,5,1);
            Ortel1 = POSfitter_vs(CRstack(:,:,1),130,TH,9,1);
            Ortel2 = POSfitter_vs(CRstack(:,:,2),130,TH,9,1);
            Ortel3 = POSfitter_vs(CRstack(:,:,3),130,TH,9,1);


            % use sum of intensity of central 3x3 pixels in PSF as quality
            % metric
            Imid_Zm = [];
            Imid_Zm(1) = mean(Ortel1(:,11));
            Imid_Zm(2) = mean(Ortel2(:,11));
            Imid_Zm(3) = mean(Ortel3(:,11));

            % use sigma of PSF as quality metric
            sig_Zm = [];
            sig_Zm(1) = mean(mean(Ortel1(:,6:7)'));
            sig_Zm(2) = mean(mean(Ortel2(:,6:7)'));
            sig_Zm(3) = mean(mean(Ortel3(:,6:7)'));

            % use symmetry of PSF as quality metric
            asym_Zm = mean(Ortes1(:,2:3) - Ortel1(:,2:3),1);
            asym_Zm = [asym_Zm; mean(Ortes2(:,2:3) - Ortel2(:,2:3),1)];
            asym_Zm = [asym_Zm; mean(Ortes3(:,2:3) - Ortel3(:,2:3),1)];

            % decide on qality metric for fit
            % Qm_Zm = (Imid_Zm)./(Ws.*sig_Zm)./(Wa+sqrt((asym_Zm(1)).^2+(asym_Zm(2)).^2));
            Qm_Zm = Wi.*Imid_Zm./10000 + (Ws*1000)./sig_Zm + (Wa*10)./(al+sqrt((asym_Zm(1)).^2+(asym_Zm(2)).^2));
            % Qm_Zm = (Imid_Zm)./sig_Zm;
            % Qm = sqrt(Imid)./sig;

            % do least-squares fit with quadratic function
            x = [-a(Zm)+p(Zm) p(Zm) a(Zm)+p(Zm)];
            pfit = polyfit(x,Qm_Zm,2);
            xfit = -a(Zm)+p(Zm):0.01:a(Zm)+p(Zm);
            yfit = pfit(1).*xfit.^2+pfit(2).*xfit+pfit(3);
            Zmfit = xfit(find(yfit==max(yfit)));
            if Zmfit<a(Zm)+p(Zm)&Zmfit>-a(Zm)+p(Zm);
                p_i(Zm) = Zmfit;
            else
                display(['Problem with mode ' num2str(Zm) ' in iteration ' num2str(i)])
                % Imid
                % sig
                % Qm
                % Zmfit
            end 

            barh(barhandle,p_i(1:21));
            xlim(barhandle,[-1 1])
            ylim(barhandle,[3.5 21+0.5])
            xlabel(barhandle,'amplitude (µm RMS)');
            ylabel(barhandle,'Zernike mode');
            yticks(barhandle,4:21)
            yticklabels(barhandle,Zyticklabels(4:21));

            if summode==1
                Zsum_i = sum(Z2C_max.*p_i(1:21),1);
                dm.Send(DMscaling(Zsum_i));
                pause(0.02)
                img = MM_getimg(mmc);
                Ortes = POSfitter_vs(img,130,TH,5,1);
                Ortel = POSfitter_vs(img,130,TH,9,1);
                Imid_i = mean(Ortel(:,11));
                sig_i = mean(mean(Ortel(:,6:7)'));
                asym_i = mean(Ortes(:,2:3) - Ortel(:,2:3),1);
                % Qm_i = (Imid_i)./(Ws.*sig_i)./(Wa+sqrt((asym_i(1)).^2+(asym_i(2)).^2));
                Qm_i = Wi.*Imid_i./10000 + (Ws*1000)./sig_i + (Wa*10)./(al+sqrt((asym_i(1)).^2+(asym_i(2)).^2));
                % Qm_i = (Imid_i)/sig_i;
                % Qm_i = sqrt(Imid_i)/sig_i;
                if Qm_i > Qm_best
                    p = p_i
                    Qm_best = Qm_i
                    Zsum = sum(Z2C_max.*p(1:21),1);
                    % a = a.*0.5;
                else
                    p_i = p;
                    if blcorr==1
                        % p_i = p;

                        % check quality metric again with previous best setting to
                        % update Qm_best for bleaching or drift, noise and fit
                        % inaccuracies
                        dm.Send(DMscaling(Zsum));
                        pause(0.02)
                        img = MM_getimg(mmc);
                        Ortes = POSfitter_vs(img,130,TH,5,1);
                        Ortel = POSfitter_vs(img,130,TH,9,1);
                        Imid_i = mean(Ortel(:,11));
                        sig_i = mean(mean(Ortel(:,6:7)'));
                        asym_i = mean(Ortes(:,2:3) - Ortel(:,2:3),1);
                        % Qm_i = (Imid_i)./(Ws.*sig_i)./(Wa+sqrt((asym_i(1)).^2+(asym_i(2)).^2));
                        Qm_i = Wi.*Imid_i./10000 + (Ws*1000)./sig_i + (Wa*10)./(al+sqrt((asym_i(1)).^2+(asym_i(2)).^2));
                        % Qm_i = (Imid_i)/sig_i;
                        Qm_best = Qm_i
                    end
                end

                barh(barhandle,p_i(1:21));
                xlim(barhandle,[-1 1])
                ylim(barhandle,[3.5 21+0.5])
                xlabel(barhandle,'amplitude (µm RMS)');
                ylabel(barhandle,'Zernike mode');
                yticks(barhandle,4:21)
                yticklabels(barhandle,Zyticklabels(4:21));
            end
        else
        end
    end
        if summode==0
            Zsum_i = sum(Z2C_max.*p_i(1:21),1);
            dm.Send(DMscaling(Zsum_i));
            pause(0.02)
            img = MM_getimg(mmc);
            Ortes = POSfitter_vs(img,130,TH,5,1);
            Ortel = POSfitter_vs(img,130,TH,9,1);
            Imid_i = mean(Ortel(:,11));
            sig_i = mean(mean(Ortel(:,6:7)'));
            asym_i = mean(Ortes(:,2:3) - Ortel(:,2:3),1);
            % Qm_i = (Imid_i)./(Ws.*sig_i)./(Wa+sqrt((asym_i(1)).^2+(asym_i(2)).^2));
            Qm_i = Wi.*Imid_i./10000 + (Ws*1000)./sig_i + (Wa*10)./(al+sqrt((asym_i(1)).^2+(asym_i(2)).^2));
            % Qm_i = (Imid_i)/sig_i;
            % Qm_i = sqrt(Imid_i)/sig_i;
            if Qm_i > Qm_best
                p = p_i;
                Qm_best = Qm_i;
                % a = a.*0.5;
            else
                p_i = p;
                % p = p_s;
            end
        end
        pall = [pall, p];
    
    Imid = [Imid, Imid_i];
    sig = [sig, sig_i];
    asym = [asym; asym_i];
    Qm = [Qm, Qm_i];
    ibest = find(Qm==max(Qm));
    pbest = pall(:,ibest);
    iteration = [iteration, i];
    barh(barhandle,p(1:21));
    xlim(barhandle,[-1 1])
    ylim(barhandle,[3.5 21+0.5])
    xlabel(barhandle,'amplitude (µm RMS)');
    ylabel(barhandle,'Zernike mode');
    yticks(barhandle,4:21)
    yticklabels(barhandle,Zyticklabels(4:21));
    plot(figh1,iteration,Imid,'o-');
    yline(figh1,Imid_start,'k--');
    yline(figh1,max(Imid),'m--');
%     plot(figh1,iteration,Imid./Imid_start,'o-');
%     yline(figh1,Imid_start./Imid_start,'k--');
%     yline(figh1,max(Imid)./Imid_start,'m--');
    xline(figh1,ibest-1,'g--');
    xlabel(figh1,'iteration')
    ylabel(figh1,'central intensity')
    xlim(figh1,[0 iter])
    xlabel(figh2,'iteration')
    ylabel(figh2,'PSF width')
    plot(figh2,iteration,sig.*2.355,'o-');
    yline(figh2,sig_start.*2.355,'k--');
    yline(figh2,min(sig).*2.355,'m--');
%      plot(figh2,iteration,sig./sig_start,'o-');
%     yline(figh2,sig_start./sig_start,'k--');
%     yline(figh2,min(sig)./sig_start,'m--');
    xline(figh2,ibest-1,'g--');
    xlim(figh2,[0 iter])
    xlabel(figh2,'iteration')
    ylabel(figh2,'PSF width')
    plot(figh3,iteration,asym(:,2),'o-');
    hold(figh3,'on')
    plot(figh3,iteration,asym(:,1),'o-');
    plot(figh3,iteration,sqrt(asym(:,1).^2+asym(:,2).^2),'o-');
    hold(figh3,'off')
    xlim(figh3,[0 iter])
    xlabel(figh3,'iteration')
    ylabel(figh3,'PSF asymmetry (nm)')
    yline(figh3,0,'k--');
    xline(figh3,ibest-1,'g--');
%     yline(figh3,asym_start,'k--');
%     atot = sqrt(asym(1,:).^2 + asym(2,:).^2);
%     amin = find(atot==min(atot));
%     amin = mean(asym(amin,:));
%     yline(figh3,amin,'m--');
    legend(figh3,'x','y','xy','Location','NorthWest')
%     plot(figh4,iteration,Qm./Qm_start,'o-');
%     yline(figh4,Qm_start./Qm_start,'k--');
%     yline(figh4,max(Qm./Qm_start),'g--');
    plot(figh4,iteration,Qm,'o-');
    yline(figh4,Qm_start,'k--');
    yline(figh4,max(Qm),'g--');
    xline(figh4,ibest-1,'g--');
    xlim(figh4,[0 iter])
    xlabel(figh4,'iteration')
    ylabel(figh4,'overall quality')
    
  
    a = a.*0.95;
end

barh(barhandle,pbest(1:21));
xlim(barhandle,[-1 1])
ylim(barhandle,[3.5 21+0.5])
xlabel(barhandle,'amplitude (µm RMS)');
ylabel(barhandle,'Zernike mode');
yticks(barhandle,4:21)
yticklabels(barhandle,Zyticklabels(4:21));

saveas(fighm,['AO_Data_' num2str(pri_only) '_pri_Zmodes_only_' num2str(prifirst) '_pri_first_' num2str(iter) 'iter_' num2str(summode) '_SumMode_' savename '_Zmodes'],'fig')
saveas(fighm,['AO_Data_' num2str(pri_only) '_pri_Zmodes_only_' num2str(prifirst) '_pri_first_' num2str(iter) 'iter_' num2str(summode) '_SumMode_' savename '_Zmodes'],'png')


Zsum = sum(Z2C_max.*pbest(1:21),1);
dm.Send(DMscaling(Zsum));

% img = MM_getimg(mmc);
% pause(0.01)
img = MM_getimg(mmc);
img_end = img;
figure('Position',[130 130 1600 1000])
subplot(2,2,1)
imagesc(img_start, [min([min(min(img_start)) min(min(img_end))])...
    max([max(max(img_start)) max(max(img_end))])]);
axis equal
colorbar
subplot(2,2,2)
imagesc(img_end, [min([min(min(img_start)) min(min(img_end))])...
    max([max(max(img_start)) max(max(img_end))])]);
axis equal
colorbar
subplot(2,2,3)
imagesc(img_start,[min([min(min(img_start)) min(min(img_end))])...
    max([max(max(img_start)) max(max(img_end))]).*0.1]);
axis equal
colorbar
subplot(2,2,4)
imagesc(img_end,[min([min(min(img_start)) min(min(img_end))])...
    max([max(max(img_start)) max(max(img_end))]).*0.1]);
axis equal
colorbar

saveas(gcf,['AO_Data_' num2str(pri_only) '_pri_Zmodes_only_' num2str(prifirst) '_pri_first_' num2str(iter) 'iter_' num2str(summode) '_SumMode_' savename '_img_start_end'],'fig')
saveas(gcf,['AO_Data_' num2str(pri_only) '_pri_Zmodes_only_' num2str(prifirst) '_pri_first_' num2str(iter) 'iter_' num2str(summode) '_SumMode_' savename '_img_start_end'],'png')


save(['AO_Data_' num2str(pri_only) '_pri_Zmodes_only_' num2str(prifirst) '_pri_first_' num2str(iter) 'iter_' num2str(summode) '_SumMode_' savename '.mat'],...
    'pall','pbest','Zsum','img_start','img_end','start_a','Ws','Wa','Wi','blcorr');

        
    



