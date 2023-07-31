import os
from datetime import datetime
from typing import Optional

import SimpleITK as sitk
import numpy as np
import click
from pydicom import uid


def gaussian(image: sitk.Image, sx: float, sy: float, sz: float) -> sitk.Image:
    # Calculate the initial resolution.
    spacing = np.array(image.GetSpacing())
    initial_sigma = spacing / (2 * np.sqrt(np.pi))

    # Calculate the filter size necessary to increase the PSF width to the given resolution. Clip at a minimum of 0
    # since we can't 'deblur' the image.
    sigma = np.sqrt((np.array((sx, sy, sz)) ** 2 - initial_sigma ** 2).clip(min=0)) / spacing

    if np.all(sigma < 0.2):
        # Kernel sizes this small have no effect: exp(-0.5 / 0.2 ** 2) < 1e-5
        return image

    gaussian_filter = sitk.DiscreteGaussianImageFilter()
    gaussian_filter.UseImageSpacingOff()  # We calculated the kernel size in pixels.
    gaussian_filter.SetVariance(sigma.astype(np.double))
    return gaussian_filter.Execute(image)


def resample(image: sitk.Image, dx: float, dy: float, dz: float) -> sitk.Image:
    # Calculate the new image size.
    size = np.array(image.GetSize())
    spacing = np.array(image.GetSpacing())
    nx, ny, nz = map(int, np.floor((size - 1) * spacing / (dx, dy, dz) + 1))  # Use floor to prevent extrapolation.

    # Resample the image. Don't change the origin and orientation.
    resample_filter = sitk.ResampleImageFilter()  # Uses linear interpolation by default.
    resample_filter.SetOutputSpacing((dx, dy, dz))
    resample_filter.SetSize((nx, ny, nz))
    resample_filter.SetOutputDirection(image.GetDirection())
    resample_filter.SetOutputOrigin(image.GetOrigin())

    return resample_filter.Execute(image)


@click.command()
@click.argument('source', type=click.Path(exists=True, file_okay=False))
@click.argument('destination', type=click.Path(dir_okay=True))
@click.option('--dx', default=None, type=float, help='Pixel spacing in x-direction')
@click.option('--dy', default=None, type=float, help='Pixel spacing in y-direction')
@click.option('--dz', default=None, type=float, help='Pixel spacing in z-direction')
def resample_dicom(source: str, destination: str, dx: Optional[float], dy: Optional[float], dz: Optional[float]):
    print("Reading DICOM series from:", source)
    reader = sitk.ImageSeriesReader()
    dicom_names = reader.GetGDCMSeriesFileNames(source)
    reader.SetFileNames(dicom_names)
    reader.MetaDataDictionaryArrayUpdateOn()  # Read and keep all meta-data.
    reader.LoadPrivateTagsOn()  # Read and keep private (odd) tags.
    image = reader.Execute()

    spacing = image.GetSpacing()
    dx = dx if dx else spacing[0]
    dy = dy if dy else spacing[1]
    dz = dz if dz else spacing[2]

    print("Smoothing volume...")
    image = gaussian(image, dx / 3.0, dy / 3.0, dz / 3.0)  # Set standard deviation to a third of the new spacing.

    print(f"Resampling volume to: {dx:.3f} × {dx:.3f} × {dz:.3f}".format(dx=dx, dy=dy, dz=dz))
    image = resample(image, dx, dy, dz)

    # Get all DICOM tags
    m = image.GetDirection()
    all_tags = {tag: reader.GetMetaData(0, tag) for tag in reader.GetMetaDataKeys(0)}
    now = datetime.now()
    all_tags["0008|0021"], datetime.strftime(now, "%Y%m%d"),  # Series Date
    all_tags["0008|0031"], datetime.strftime(now, "%H%M%S"),  # Series Time
    all_tags["0008|0008"] = all_tags.get("0008|0008", "ORIGINAL\\PRIMARY").replace("ORIGINAL", "DERIVED")  # Image Type
    all_tags["0020|000e"] = uid.generate_uid()  # Series Instance UID
    all_tags["0020|0037"] = "\\".join(map(str, (m[0], m[3], m[6], m[1], m[4], m[7])))  # Image Orientation (Patient)

    # Remove optional slice specific tags.
    del all_tags["0020|1041"]  # Slice Location

    # Set new Slice Thickness.
    if "0018|0050" in all_tags:
        all_tags["0018|0050"] = f"{dz:.3f}".format(dz=dz)

    writer = sitk.ImageFileWriter()
    writer.KeepOriginalImageUIDOn()  # Use the study/series/frame of reference information

    print("Writing DICOM series to:", destination)
    for i in range(image.GetDepth()):
        image_slice = image[:, :, i]

        # Set DICOM tags.
        for tag, value in all_tags.items():
            image_slice.SetMetaData(tag, value)

        # Set slice specific tags.
        now = datetime.now()
        instance_uid = uid.generate_uid()
        pos = image.TransformIndexToPhysicalPoint((0, 0, i))
        image_slice.SetMetaData("0008|0012", datetime.strftime(now, "%Y%m%d"))  # Instance Creation Date
        image_slice.SetMetaData("0008|0013", datetime.strftime(now, "%H%M%S"))  # Instance Creation Time
        image_slice.SetMetaData("0008|0018", instance_uid)  # SOP Instance UID
        image_slice.SetMetaData("0020|0032", "\\".join(map(str, pos)))  # Image Position (Patient)
        image_slice.SetMetaData("0020|0013", str(i))  # Instance Number

        # Write to the output directory and add the extension dcm, to force writing
        # in DICOM format.
        if not os.path.exists(destination):
            os.makedirs(destination)
        writer.SetFileName(os.path.join(destination, instance_uid + ".dcm"))
        writer.Execute(image_slice)


if __name__ == '__main__':
    resample_dicom()
