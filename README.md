how to run server:

1. set up venv, so that:
python -m venv nazwa_venv


folder:
    -ur_venv
    -cos.py


2. activate venv
nazwa_venv\Scripts\activate

3. download stuff
pip install fastapi uvicorn

4. run
uvicorn cos:app --reload


5. open this in browser
http://127.0.0.1:8000/docs


6. close
ctrl + c

7. deactivate venv
deactivate
