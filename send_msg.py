import sys
import requests
import os
import configparser

def load_config():
    config_path = os.path.expanduser("~/telegram_config.ini")  # 配置文件路径
    config = configparser.ConfigParser()
    config.read(config_path)
    return config['Telegram']['bot_token'], config['Telegram']['chat_id']

def send_message_to_telegram(message):
    bot_token, chat_id = load_config()

    url = f"https://api.telegram.org/bot{bot_token}/sendMessage"
    params = {
        "chat_id": chat_id,
        "text": message
    }

    response = requests.post(url, json=params)
    if response.status_code == 200:
        print("Message sent successfully.")
    else:
        print("Failed to send message.")

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python send_msg.py <message>")
        sys.exit(1)

    message = sys.argv[1]
    send_message_to_telegram(message)