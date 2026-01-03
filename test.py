from transformers import T5Tokenizer, T5ForConditionalGeneration
import torch

# Load model and tokenizer from local directory
MODEL_PATH = "models/flan-t5-base"  # Make sure this path exists and contains model files
tokenizer = T5Tokenizer.from_pretrained(MODEL_PATH)
model = T5ForConditionalGeneration.from_pretrained(MODEL_PATH)

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
        max_length=500,
        temperature=0.9,
        top_k=50,
        top_p=0.95,
        do_sample=True,
        num_return_sequences=1,
        repetition_penalty=1.3
    )
    explanation = tokenizer.decode(output[0], skip_special_tokens=True)
    return explanation


# 🔧 Take user input correctly
molecule = input("Enter molecule name (e.g., ethanol): ").strip().lower()
print("\n🔬 Explanation:")
print(explain_molecule(molecule))
