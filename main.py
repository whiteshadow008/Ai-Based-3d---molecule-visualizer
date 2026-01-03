import os
import time
import traceback
from flask import Flask, render_template, request, jsonify, send_from_directory

from generate_xyz import generate_xyz
from blender_script_gen import generate_blender_py
from export_glb_blender import run_blender_script

# LLM (Offline)
from transformers import T5Tokenizer, T5ForConditionalGeneration
import torch

# Load model and tokenizer from local directory
MODEL_PATH = "models/flan-t5-base"
tokenizer = T5Tokenizer.from_pretrained(MODEL_PATH)
model = T5ForConditionalGeneration.from_pretrained(MODEL_PATH)
from transformers import T5Tokenizer, T5ForConditionalGeneration # type: ignore

tokenizer = T5Tokenizer.from_pretrained("models/flan-t5-base")
model = T5ForConditionalGeneration.from_pretrained("models/flan-t5-base")

app = Flask(__name__)

# Molecule SMILES dictionary
molecule_db = {
    "water":"H2O",
    "methane": "C", "ethane": "CC", "propane": "CCC", "butane": "CCCC", "pentane": "CCCCC",
    "hexane": "CCCCCC", "ethene": "C=C", "ethyne": "C#C", "benzene": "c1ccccc1", "methanol": "CO",
    "ethanol": "CCO", "propanol": "CCCO", "butanol": "CCCCO", "isopropanol": "CC(O)C",
    "glycerol": "C(C(CO)O)O", "formic acid": "C(=O)O", "acetic acid": "CC(=O)O", "propanoic acid": "CCC(=O)O",
    "butanoic acid": "CCCC(=O)O", "benzoic acid": "c1ccc(cc1)C(=O)O", "methylamine": "CN",
    "ethylamine": "CCN", "aniline": "c1ccccc1N", "formaldehyde": "C=O", "acetaldehyde": "CC=O",
    "acetone": "CC(=O)C", "glucose": "OC[C@H]1O[C@@H](O)[C@H](O)[C@@H](O)[C@H]1O",
    "fructose": "C(C1C(C(C(C(O1)(CO))O)O)O)=O", "sucrose": "OCC1OC(O)C(O)C(O)C1OC2(C(C(C(C(O2)CO)O)O)O)",
    "glycine": "NCC(=O)O", "alanine": "CC(C(=O)O)N", "valine": "CC(C)C(C(=O)O)N", "leucine": "CC(C)CC(C(=O)O)N",
    "phenylalanine": "C1=CC=C(C=C1)CC(C(=O)O)N", "tryptophan": "C1=CC=C2C(=C1)C(=CN2)CC(C(=O)O)N",
    "ammonia": "N", "carbon dioxide": "O=C=O", "nitric acid": "O=N(=O)O", "sulfuric acid": "OS(=O)(=O)O",
    "phosphoric acid": "OP(=O)(O)O", "hydrogen peroxide": "OO", "ozone": "O=O[O]",
    "adenine": "C1=NC2=C(N1)N=CN2", "guanine": "C1=NC2=C(N1)C(=O)N=CN2", "cytosine": "C1=CN=CN1",
    "uracil": "C1=CC(=O)NC(=O)N1", "thymine": "CC1=CN(C(=O)NC1=O)C", 
    "cholesterol": "C[C@H](CCC(=O)O)C1CCC2C3CCC4=CC(=O)CCC4(C)C3CCC12C",
    "caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C", "nicotine": "CN1CCCC1C2=CN=CC=C2",
    "oxygen": "O=O", "hydrogen": "[H][H]", "nitrogen": "N#N", "chlorine": "ClCl", "fluorine": "F[F]",
    "aspirin": "CC(=O)Oc1ccccc1C(=O)O", "paracetamol": "CC(=O)NC1=CC=C(C=C1)O",
    "acetylsalicylic acid": "CC(=O)Oc1ccccc1C(=O)O"
}

# Molecule explanation generator using LLM
def explain_molecule(molecule: str) -> str:
    prompt = (
        f"You are a chemistry expert. Write a detailed explanation about the molecule '{molecule}'that a school student can easily understand, including:\n"
        f"- Its molecular structure and chemical formula\n"
        f"- Functional groups present\n"
        f"- Physical and chemical properties\n"
        f"- Common uses in real life or industry\n"
        f"- Any biological or environmental relevance\n\n"
        f"Start your explanation now:"
    )

    input_ids = tokenizer(prompt, return_tensors="pt").input_ids
    output = model.generate(
        input_ids,
        max_length=400,
        temperature=0.9,
        top_k=50,
        top_p=0.95,
        do_sample=True,
        num_return_sequences=1,
        repetition_penalty=1.3
    )
    explanation = tokenizer.decode(output[0], skip_special_tokens=True)
    return explanation

@app.route('/', methods=['GET', 'POST'])
def index():
    if request.method == 'POST':
        try:
            name = request.form['molecule'].strip().lower()
            if name in molecule_db:
                smiles = molecule_db[name]
                os.makedirs('outputs', exist_ok=True)

                xyz_path = generate_xyz(smiles, name)
                blender_script = generate_blender_py(xyz_path, name)
                run_blender_script(blender_script)

                glb_path = f"outputs/{name}.glb"
                for _ in range(40):  # wait up to 20 seconds
                    if os.path.exists(glb_path):
                        print(f"✅ File created: {glb_path}")
                        return jsonify(success=True, viewer_url=f"/viewer/{name}")
                    time.sleep(0.5)

                print(f"❌ GLB file not created: {glb_path}")
                return jsonify(success=False, error="Failed to generate 3D model.")
            else:
                return jsonify(success=False, error="Invalid molecule name.")
        except Exception as e:
            print("❌ Error:")
            traceback.print_exc()
            return jsonify(success=False, error=str(e))
    return render_template('index.html')

@app.route('/viewer/<molecule>')
def viewer(molecule):
    glb_file = f"outputs/{molecule}.glb"
    if not os.path.exists(glb_file):
        return f"File not found: {glb_file}", 404

    explanation = explain_molecule(molecule)
    return render_template("viewer.html", molecule=molecule, explanation=explanation)

@app.route('/outputs/<path:filename>')
def serve_outputs(filename):
    return send_from_directory("outputs", filename)

if __name__ == '__main__':
    app.run(debug=True, use_reloader=False)