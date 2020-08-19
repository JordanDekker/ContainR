FROM python:3.6

# Prepare app for deployment
WORKDIR /containr
COPY ./LICENSE ./requirements.txt ./webapp.py ./config.py ./install_blast.sh ./
ENV PATH "$PATH:/containr/ncbi-blast-2.9.0+/bin"
RUN sh install_blast.sh
RUN pip install -r requirements.txt
COPY ./app ./app

EXPOSE 5000

CMD ["gunicorn", "-b", ":5000", "--access-logfile", "-", "--error-logfile", "-", "--workers=13", "--threads=6", "webapp:app"]

# CMD echo "Container is now running, press Ctrl-C to stop the container." & tail -f /dev/null