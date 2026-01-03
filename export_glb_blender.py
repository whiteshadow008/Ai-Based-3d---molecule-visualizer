import subprocess

def run_blender_script(script_path):
    blender_path = r"C:\Program Files\Blender Foundation\Blender 4.4\blender.exe"  # ✅ Update path if needed

    cmd = [
        blender_path,
        "--background",
        "--python",
        script_path
    ]

    try:
        result = subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        print("✅ Blender Output:\n", result.stdout)
        if result.stderr:
            print("⚠️ Blender Errors:\n", result.stderr)
    except subprocess.CalledProcessError as e:
        print("❌ Blender Failed:")
        print("STDOUT:\n", e.stdout)
        print("STDERR:\n", e.stderr)
        raise  # Let it bubble up to the main app for full traceback
