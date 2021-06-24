function fov_guider

directory = '/home/jaehan/Desktop/DJI_202104281514_048';
fileList = dir(directory);
fileList = fileList(3:end);
fileList.name;

listSize = size(fileList,1)

while 1
    imageNum = input('Enter Image Number');
    
    zoom = 5;
    targetZoom = 20;
    ratio = targetZoom / zoom;

    image = imread(directory + "/" + fileList(imageNum).name);
    hpixel = size(image,2);
    vpixel = size(image,1);
    pixelCenter = [hpixel/2,vpixel/2];
    guideH = hpixel / ratio / 2;
    guideV = vpixel / ratio / 2;

    figure(10)
    imshow(image)
    hold on
%     plot([pixelCenter(1)-guideH pixelCenter(1)+guideH],[pixelCenter(2)-guideV pixelCenter(2)-guideV],'r','LineWidth',4)
%     plot([pixelCenter(1)-guideH pixelCenter(1)+guideH],[pixelCenter(2)+guideV pixelCenter(2)+guideV],'r','LineWidth',4)
%     plot([pixelCenter(1)+guideH pixelCenter(1)+guideH],[pixelCenter(2)-guideV pixelCenter(2)+guideV],'r','LineWidth',4)
%     plot([pixelCenter(1)-guideH pixelCenter(1)-guideH],[pixelCenter(2)-guideV pixelCenter(2)+guideV],'r','LineWidth',4)
    plot([pixelCenter(1) pixelCenter(1)],[0 vpixel],'r','LineWidth',4)
end

end