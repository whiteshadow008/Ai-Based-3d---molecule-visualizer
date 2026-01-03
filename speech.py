from gtts import gTTS
import os

# Heart explanation text
heart_text = """The human skeletal system forms the structural framework of the body.
It consists of 206 bones in an adult human.
Bones provide support and give shape to the body.
They protect delicate organs like the brain, heart, and lungs.
The skeletal system also allows movement when muscles pull on bones.
Bones store minerals such as calcium and phosphorus.
They also produce blood cells in the bone marrow.
The skeleton is divided into two main parts: axial and appendicular.
The axial skeleton includes the skull, vertebral column, and rib cage.
It protects the brain, spinal cord, and thoracic organs.
The appendicular skeleton includes arms, legs, shoulder, and hip bones.
It helps in locomotion and manipulation of the environment.
Joints connect bones and allow movement.
There are different types of joints like hinge, ball-and-socket, and pivot.
Ligaments hold bones together, while tendons connect muscles to bones.
The skull has 22 bones, including the jaw (mandible).
The vertebral column has 33 vertebrae arranged in regions.
The rib cage consists of 12 pairs of ribs protecting the heart and lungs.
Bone health depends on calcium, vitamin D, and regular exercise.
Overall, the skeletal system is vital for support, protection, and movement.
"""

def text_to_audio(text, filename="skeleton.mp3", lang="en"):
    try:
        # Get user Downloads folder path
        downloads_path = os.path.join(os.path.expanduser("~"), "Downloads")
        
        # Ensure folder exists
        if not os.path.exists(downloads_path):
            os.makedirs(downloads_path)
        
        # Full file path inside Downloads
        file_path = os.path.join(downloads_path, filename)
        
        # Convert text to speech and save
        tts = gTTS(text=text, lang=lang, slow=False)
        tts.save(file_path)
        
        print(f"✅ Audio saved in Downloads: {file_path}")
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    text_to_audio(heart_text)
