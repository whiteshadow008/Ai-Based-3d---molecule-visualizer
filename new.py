from transformers import T5Tokenizer, T5ForConditionalGeneration

# Use a pipeline as a high-level helper
from transformers import pipeline

pipe = pipeline("translation", model="google-t5/t5-small")

print("T5 Interactive Mode (type 'quit' to exit)")

while True:
    text = input("Enter input: ")
    if text.lower() == "quit":
        break
    
    inputs = tokenizer(text, return_tensors="pt")
    output_ids = model.generate(inputs["input_ids"], max_length=50)
    output = tokenizer.decode(output_ids[0], skip_special_tokens=True)
    print("Model Output:", output)
