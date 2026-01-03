from transformers import T5Tokenizer, T5ForConditionalGeneration

model_name = "google/flan-t5-base"
save_path = "models/flan-t5-base"  # Local folder

# Download and save tokenizer
tokenizer = T5Tokenizer.from_pretrained(model_name)
tokenizer.save_pretrained(save_path)

# Download and save model
model = T5ForConditionalGeneration.from_pretrained(model_name)
model.save_pretrained(save_path)

print("✅ Model and tokenizer saved to:", save_path)
