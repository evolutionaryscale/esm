import textwrap


def wrapped_print(text, width=70):
    text = str(text)
    wrapped_text = textwrap.fill(text, width=width)
    print(wrapped_text)
