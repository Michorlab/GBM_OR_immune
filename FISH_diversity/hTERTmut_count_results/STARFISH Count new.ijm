v = getInfo("image.filename"); //this first bit takes the name of the file and defines the variables needed for renaming all the new files//
l = lengthOf(v);
v1 = substring(v, 0, l-4);
v2 = substring(v, 0, l-5);
v3 = substring(v, 0, l-8);
p = getDirectory("image");
run("Duplicate...", v);
selectWindow(v1+"-1.tif");
run("8-bit");
setAutoThreshold("Default dark");
waitForUser("Adjust brightness and threshold. Click OK when ready");
run("Convert to Mask");
waitForUser("Fill Holes? Click OK when ready?");
run("Watershed");
run("Analyze Particles...", "size=200-Infinity pixel circularity=0.00-1.00 show=Outlines clear include add");
selectWindow("Drawing of "+v1+"-1.tif");
save(p+"Drawing of "+v1+".tif");
close();
roiManager("Show All with labels");


selectWindow(v2+"3.tif"); //needs tresholding
run("8-bit");
waitForUser("Adjust brightness and threshold. Click OK when ready");
run("Clear Results");
//run("Invert");
n = roiManager("count");
for (i=0; i<n; i++) {
     roiManager("select", i);
run("Analyze Particles...", "size=4-200 pixel circularity=0.00-1.00 show=Nothing display");
roiManager("deselect");
print(i+1, ",", nResults);
run("Clear Results");
}

selectWindow("Log");
save(p+"Log_"+v3+"tert.txt");
//saveAs("Text", "/Users/katherinemurphy/Desktop/STARFISH/Log_"+v2+"2.txt");
//saveAs("Text", "/Users/katherinemurphy/Desktop/STARFISH/TDM1/log/Log_"+v2+"2.txt");

print("\\Clear");


roiManager("save", p+"Roi_"+v1+".zip");

selectWindow(v1+"-1.tif");
run("Set Measurements...", "area centroid limit display redirect=None decimal=0");
n = roiManager("count");
for(i=0; i<n; i++) {
	roiManager("select", i);
	run("Measure");
	x = getResult('X');
	y = getResult('Y');
	a = getResult("Area");
getVoxelSize(width, height, depth, unit);
x1 = x/width;
y1 = y/height;
	print(i+1, ",", d2s(x1,0), ",", d2s(y1,0), ",", a); //d2s function is for decimal places d2s(pi,2) => 3.14
run("Clear Results");
}

selectWindow("Log");
save(p+"Log_"+v3+"_xyarea.txt");
//saveAs("Text", "/Users/katherinemurphy/Desktop/STARFISH/Log_"+v1+"_xy.txt");
//saveAs("Text", "/Users/katherinemurphy/Desktop/STARFISH/TDM1/log/Log_"+v1+"_xy.txt");

roiManager("Deselect");
roiManager("Delete");
run("Close");