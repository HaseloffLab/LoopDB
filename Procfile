#web: gunicorn --worker-class eventlet -w 1 q server:app
#web: gunicorn server:app
#web: uwsgi --http :5000 --gevent 1000 --http-websockets --master --wsgi-file server.py --callable app
web: gunicorn -k geventwebsocket.gunicorn.workers.GeventWebSocketWorker -w 1 server:app
