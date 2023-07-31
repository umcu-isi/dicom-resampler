# DICOM Resampler
Resamples 3D DICOM series. The tool can be use for both upsampling and downsampling 3D volumes to a specified pixel spacing. Anisotropic Gaussian smoothing is applied when necessary.

## Requirements
```shell
pip install SimpleITK numpy click pydicom
```

## Example
Create 3mm thick axial slices. Read from /input/files and write to /output/resampled: 
```shell
python resample.py /input/files /output/resampled --dz 3
```

Create 1mm isotropic voxels. Read from /input/files and write to /output/resampled: 
```shell
python resample.py /input/files /output/resampled --dx 1 --dy 1 --dz 1
```
