function  Subtract_then_fit(mov_fname,Mol_off_frames_fname,guessfname,MLE_fit,edgedist,stdtol,maxerr,do_avgsub,which_gaussian)
%% Subtract_mol_off_frames
% subtracts the average (or median) intensity of off frames for each guess
% stored in Mol_off_frames_fname.
%
% If you just want to do fitting, and not do background subtraction, set
% Mol_off_frames_fname = 'nobgsub'. The program will take care of
% everything else.

%%%% Inputs %%%%
% mov_fname the filename of the tiff stack movie

% Mol_off_frames_fname is the filename .mat file output from the function
% Mol_off_frames. If not doing background subtraction, set this to
% 'nobgsub'

% guessfname is the filename for the guesses .mat file

% MLE_fit  a Boolean determining whether or not MLE fitting is used. Set to
% 1 to use MLE and to 0 to use least squares. Default is 0. Note that MLE
% is quite slow, and so its not recommended for a large number of guesses

% edgedist is the distance in pixels from the edge of the frame to ignore.
% default is 10

% stdtol is tolerance on fit Gaussian STD, to leae filtering options for
% later, default value is 1.5

% maxerr is the maximum error of the fit for MLE fit, using variance default
% 0.1 (can't be above this) for LSQR fit, using the 95% confidence interval
% on the position, default max is 2

% do_avgsub is a Boolean determining whether or not to subtract the mean of
% the off frames. Set to 1 to subtract the mean and to 0 to subtract the
% median. Default is 1.

% which_gaussian determines what functional form of Gaussian function the
% molecules will be fit to if using least-squares fitting (MLE fitting only
% fits symmetric Gaussian). Set to 1 to use a symmetric Gaussian. Set to 2
% to use an asymmetric Gaussian (with axes along the row and column
% dimension). Set to 3 to use a freely rotating asymmetric Gaussian.
% Default is 1.

%%%% Output %%%%
% a .mat file, importantly containing the fits structure that has fields
%frame number of the fit:
% fits.frame
%row coordinate of the fit:
% fits.row
%column coordinate of the fit:
% fits.col
%standard deviation in the row dimension of the Gaussian fit (if using a
%symmetric Gaussian this will be the same as the other width):
% fits.widthr
%standard deviation in the column dimension of the Gaussian fit (if using a
%symmetric Gaussian this will be the same as the other width):
% fits.widthc
%angle of asymmetric Gaussian fit:
% fits.ang
%offset of Gaussian fit:
% fits.offset
%amplitude of Gaussian fit:
% fits.amp
%error on fit (for MLE fitting, this is the variance, for least squares
%fitting, this is the mean 95% confidence interval on the position):
% fits.err
%sum of pixels in ROI around guess:
% fits.sum
%goodfit boolean:
% fits.goodfit

%%%% Dependencies %%%%
% TIFFStack
% MLEwG (for MLE fitting)
% gaussfit (for least squares fitting)

%     Copyright (C) 2017  Benjamin P Isaacoff
%
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%


if nargin<4;MLE_fit=0;end
if nargin<5;edgedist=10;end
if nargin<6;stdtol=1.5;end
if nargin<7;
    if MLE_fit
        maxerr=0.1;
    else
        maxerr=2;
    end
end
if nargin<8;do_avgsub=1;end
if nargin<9;which_gaussian=1;end

tic;%for measuring the time to run the entire program
%% Import the data

%create A `TIFFStack` object  which behaves like a read-only memory
%mapped TIFF file
tfstk=TIFFStack(mov_fname);
movsz=size(tfstk);%the size of the movie
[pathstr,fname,~] = fileparts(mov_fname);

% load the guesses
load(guessfname,'guesses','dfrlmsz');

%check that the bgsub is actually happening
bgsub=1;
if strcmp(Mol_off_frames_fname,'nobgsub');bgsub=0;end

if bgsub
    % load off frames list and some parameters
    load(Mol_off_frames_fname,'off_frames','moloffwin')
else
    %set moloffwin to a fifth of the frames, just for memory purposes
    moloffwin=ceil((movsz(3)/5)/2)*2;
    %create filled off_frames cell for simplicity, this isn't used for
    %anything other than not being empty
    off_frames=cell([size(guesses,1),1]);
    off_frames(:)={'foobar'};
end

if bgsub
    %check number of fits vs length of off frames
    if size(guesses,1)~=numel(off_frames);error('Unequal number of fits and number of off frames lists');end
end

% import the first moloffwin+1 frames
curframes=1:(moloffwin+1);
mov=double(tfstk(:,:,curframes));

%% The Averaging and Subtraction

%the conversion between dfrlmsz and the STD of the Gaussian, reccomended
%using the full width at 20% max given by (2*sqrt(2*log(5)))
dfD2std=(2*sqrt(2*log(5)));
%the guessed std
gesss=dfrlmsz/dfD2std;

%initializing the fits structure
fits.frame=NaN(size(guesses,1),1);%frame numbers
fits.row=NaN(size(guesses,1),1);%row coordinate of the fit
fits.col=NaN(size(guesses,1),1);%column coordinate of the fit
fits.widthr=NaN(size(guesses,1),1);%standard deviation in the row dimension of the Gaussian fit
fits.widthc=NaN(size(guesses,1),1);%standard deviation in the column dimension of the Gaussian fit
fits.ang=NaN(size(guesses,1),1);%angle of asymmetric Gaussian fit
fits.offset=NaN(size(guesses,1),1);%offset
fits.amp=NaN(size(guesses,1),1);%amplitude of Gaussian fit
fits.err=NaN(size(guesses,1),1);%error on fit
fits.sum=NaN(size(guesses,1),1);%sum of pixels in ROI around guess
fits.goodfit=false(size(guesses,1),1);%goodfit boolean

%starting the waitbar
h1=waitbar(0);
set(findall(h1,'type','text'),'Interpreter','none');
waitbar(0,h1,['Fitting ',fname]);

%looping through all the guesses
for ii=1:size(guesses,1)
    try; waitbar(ii/size(guesses,1),h1); end
    
    %current frame number
    curfrmnum=guesses(ii,1);
    %putting it into the fits structure
    fits.frame(ii)=curfrmnum;%frame numbers
    
    %determine the frame list of frames to check for the current frame
    if curfrmnum<=(moloffwin/2)%the first group of frames
        frmlst=curfrmnum+(-(curfrmnum-1):(moloffwin/2));
    elseif curfrmnum>=(movsz(3)-moloffwin/2)%the last group of frames
        frmlst=movsz(3)+(-moloffwin:0);
    else %all the frames in the middle
        frmlst=curfrmnum+((-moloffwin/2):(moloffwin/2));
    end
    
    %import appropriate movie frames
    if frmlst(end)>curframes(end)
        numnewfrmsend=frmlst(end)-curframes(end);
        numnewfrmsbeg=frmlst(1)-curframes(1);
        %new current frames list
        curframes=frmlst;
        %new movie frames
        if numnewfrmsend<(moloffwin+1)
            mov=cat(3,mov(:,:,(numnewfrmsbeg+1):(moloffwin+1)),...
                double(tfstk(:,:,curframes(end)+(-(numnewfrmsend-1):0))));
        else
            mov=double(tfstk(:,:,curframes));
            warning('There was a big jump in the frames without a fit')
        end
        %check to make sure that everything is the right size
        if curfrmnum>(moloffwin) && (length(curframes)~=(moloffwin+1) || size(mov,3)~=(moloffwin+1))
            error('Error determining the correct frames to import')
        end
    end
    
    %current molecule's position
    molr=guesses(ii,2);
    molc=guesses(ii,3);
    
    %checking that it's not outside the frame and that off_frames for this
    %guess isn't empty
    if (molc>edgedist && molc<(movsz(2)-edgedist) && molr>edgedist && molr<(movsz(1)-edgedist)) && ...
            (molc>dfrlmsz && molc<(movsz(2)-dfrlmsz) && molr>dfrlmsz && molr<(movsz(1)-dfrlmsz))&& ...
            ~isempty(off_frames{ii})
        if bgsub
            %the average (or median) frame
            if do_avgsub
                mean_mov=mean(mov(molr+(-dfrlmsz:dfrlmsz),molc+(-dfrlmsz:dfrlmsz),off_frames{ii}-frmlst(1)+1),3);
            else
                mean_mov=median(mov(molr+(-dfrlmsz:dfrlmsz),molc+(-dfrlmsz:dfrlmsz),off_frames{ii}-frmlst(1)+1),3);
            end
            %the molecule image
            molim=mov(molr+(-dfrlmsz:dfrlmsz),molc+(-dfrlmsz:dfrlmsz),curfrmnum-frmlst(1)+1);
            %the subtracted image
            data=molim-mean_mov;
        else
            data=mov(molr+(-dfrlmsz:dfrlmsz),molc+(-dfrlmsz:dfrlmsz),curfrmnum-frmlst(1)+1);
        end
        
        %%%% Fitting %%%%
        plot_on=0;%for debugging purposes only!
        %the guessed tail intensity of the gaussian
        gessb=min(data(:));
        %the guessed amplitude, using the formula in MLEwG
        gessN=range(data(:))*(4*pi*gesss^2);
        %fit guess vector
        params0=[dfrlmsz,dfrlmsz,gesss,gessb,gessN];
        
        if MLE_fit
            %fitting with MLE
            [paramsF,varianceF] = MLEwG (data,params0,1,plot_on,1);
            %shifting
            paramsF([1,2])=paramsF([1,2])+0.5;
            fit_r=paramsF(1);fit_c=paramsF(2);
            fit_sd_r=paramsF(3);fit_sd_c=paramsF(3);
            %recalculating the values based on their equations to match
            paramsF(5)=paramsF(5)/(2*pi*paramsF(3)^2);
            if paramsF(4)>=0
                paramsF(4)=sqrt(paramsF(4));
            else
                paramsF(4)=-sqrt(-paramsF(4));
            end
            fit_off=paramsF(4);
            fit_amp=paramsF(5);
            fit_ang=0;
            fit_err=varianceF;
            errbad=varianceF>maxerr;%too much error on fit?
        else
            %fitting with least squares
            [fitPars,conf95,~,~,~]=gaussFit(data,'searchBool',0,'nPixels',2*dfrlmsz+1,...
                'checkVals',0,'ffSwitch',which_gaussian);
            %converting the variables to match the output of MLEwG, and
            %arranging for each particular Gaussian fit
            fit_r=fitPars(1);fit_c=fitPars(2);
            if which_gaussian==1
                fit_sd_r=fitPars(3);fit_sd_c=fitPars(3);
                fit_off=fitPars(5);
                fit_amp=fitPars(4);                
                fit_ang=0;
            elseif which_gaussian==2
                fit_sd_r=fitPars(3);fit_sd_c=fitPars(4);
                fit_off=fitPars(6);
                fit_amp=fitPars(5);
                fit_ang=0;
            elseif   which_gaussian==3
                fit_sd_r=fitPars(4);fit_sd_c=fitPars(5);
                fit_off=fitPars(7);
                fit_amp=fitPars(6);
                fit_ang=fitPars(3);
            end
            fit_err=mean(conf95([1,2]));
            errbad=mean(conf95([1,2]))>maxerr;%too much error on fit?
        end
        %Convert back into full frame coordinates, NOTE the -1!
        act_r=fit_r-dfrlmsz-1+molr;
        act_c=fit_c-dfrlmsz-1+molc;
        %The sum(:) of the the data
        sumsum=sum(data(:));
        
        %putting the fit results into the fits structure
        fits.row(ii)=act_r;%row coordinate of the fit
        fits.col(ii)=act_c;%column coordinate of the fit
        fits.widthr(ii)=fit_sd_r;%standard deviation in the row dimension of the Gaussian fit
        fits.widthc(ii)=fit_sd_c;%standard deviation in the column dimension of the Gaussian fit
        fits.ang(ii)=fit_ang;%angle of asymmetric Gaussian fit
        fits.offset(ii)=fit_off;%offset
        fits.amp(ii)=fit_amp;%amplitude of Gaussian fit
        fits.err(ii)=fit_err;%error on fit
        fits.sum(ii)=sumsum;%sum of pixels in ROI around guess
        
        %determining if it's a goodfit or not (remember this field was
        %initialized to false)
        if (mean([fit_sd_r,fit_sd_c])<=(stdtol*params0(3)) && mean([fit_sd_r,fit_sd_c])>=(params0(3)/stdtol)) && ... %Compare width with diffraction limit
                ~errbad && ... %too much error on fit?
                fit_amp<sumsum && ... %the amplitude of the fit shouldn't be bigger than the integral
                ~any([fit_r,fit_c,fit_amp,sumsum]<0) %none of the fitted parameters should be negative, except the offset!
            
            fits.goodfit(ii)=true;%goodfit boolean
        end
        
        %plotting for debugging/tests
        if plot_on
            h12=figure(12);
            subplot(1,3,1)
            imshow(mean_mov,[])
            title('Mean BG')
            subplot(1,3,2)
            imshow(molim,[])
            title('Raw Molecule')
            subplot(1,3,3)
            imshow(data,[])
            title('BGSUB')
            annotation('textbox', [0 0.9 1 0.1], ...
                'String', ['Frame # ',num2str(curfrmnum),'   Guess # ',num2str(ii)], ...
                'EdgeColor', 'none', ...
                'HorizontalAlignment', 'center')
            
            keyboard
            close(h12)
        end
    end
end

tictoc=toc;%the time to run the entire program
%save the data
if bgsub
    fname=[pathstr,filesep,fname,'_AccBGSUB_fits.mat'];
else
    fname=[pathstr,filesep,fname,'_fits.mat'];
end
save(fname,'fits','mov_fname','Mol_off_frames_fname','guessfname',...
    'MLE_fit','stdtol','maxerr','dfrlmsz','movsz','moloffwin','tictoc',...
    'do_avgsub','which_gaussian')

try
    close(h1)
end
end

