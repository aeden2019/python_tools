import os
from PIL import Image
import sys

# Path to the folder containing your images
image_folder = '/home/jovyan/home/images/star_subtracted_tan_vel'

# Get a list of all image files in the folder
image_files = [f for f in os.listdir(image_folder) if f.endswith('.png')]

# Sort the image files if necessary
image_files.sort()

# Create a list to store image objects
images = []

# Open each image and add it to the images list
for image_file in image_files:
    image_path = os.path.join(image_folder, image_file)
    image = Image.open(image_path)
    images.append(image)
    
# Specify the output file name for the animation
#output_file = 'animation.gif'
output_file = sys.argv[1]

# Save the animation
images[0].save(
    output_file,
    save_all=True,
    append_images=images[1:],
    duration=200,  # Duration between frames in milliseconds
    loop=0,       # Loop forever
    optimize=False  # Disable optimization for better quality
)

print("Animation saved as", output_file)