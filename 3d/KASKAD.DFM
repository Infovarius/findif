object MainWindow: TMainWindow
  Left = 210
  Top = 118
  Width = 757
  Height = 504
  Caption = 'Интегрирование каскадно-сеточным методом'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  OldCreateOrder = False
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  PixelsPerInch = 96
  TextHeight = 13
  object Notebook: TTabbedNotebook
    Left = 0
    Top = 0
    Width = 757
    Height = 479
    Anchors = [akLeft, akTop, akRight, akBottom]
    TabFont.Charset = DEFAULT_CHARSET
    TabFont.Color = clBtnText
    TabFont.Height = -11
    TabFont.Name = 'MS Sans Serif'
    TabFont.Style = []
    TabOrder = 0
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Initiation'
      object GroupBox1: TGroupBox
        Left = 296
        Top = 24
        Width = 337
        Height = 313
        Caption = 'Имена файлов вывода'
        TabOrder = 0
        object Label1: TLabel
          Left = 16
          Top = 36
          Width = 89
          Height = 13
          Caption = 'Файл сообщений'
        end
        object Label2: TLabel
          Left = 16
          Top = 64
          Width = 127
          Height = 13
          Caption = 'Вывод значений энергии'
        end
        object Label3: TLabel
          Left = 16
          Top = 96
          Width = 133
          Height = 13
          Caption = 'Вывод значений скорости'
        end
        object Label4: TLabel
          Left = 16
          Top = 156
          Width = 157
          Height = 13
          Caption = 'Вывод каскадных переменных'
        end
        object Label5: TLabel
          Left = 16
          Top = 196
          Width = 70
          Height = 13
          Caption = 'Файл ошибок'
        end
        object Label6: TLabel
          Left = 16
          Top = 260
          Width = 64
          Height = 13
          Caption = 'Файл дампа'
        end
        object Label7: TLabel
          Left = 16
          Top = 124
          Width = 133
          Height = 13
          Caption = 'Вывод значений вязкости'
        end
        object EditMessage: TEdit
          Left = 152
          Top = 32
          Width = 145
          Height = 21
          TabOrder = 0
          Text = 'message.dat'
          OnKeyUp = EditMessageKeyUp
        end
        object EditEnergy: TEdit
          Left = 152
          Top = 60
          Width = 145
          Height = 21
          TabOrder = 1
          Text = 'energy.dat'
          OnKeyUp = EditEnergyKeyUp
        end
        object EditNu: TEdit
          Left = 152
          Top = 120
          Width = 145
          Height = 21
          TabOrder = 3
          Text = 'nut.dat'
          OnKeyUp = EditNuKeyUp
        end
        object EditKaskvar: TEdit
          Left = 184
          Top = 152
          Width = 113
          Height = 21
          TabOrder = 4
          Text = 'kaskvar.dat'
          OnKeyUp = EditKaskvarKeyUp
        end
        object EditError: TEdit
          Left = 152
          Top = 192
          Width = 145
          Height = 21
          TabOrder = 5
          Text = 'error.err'
          OnKeyUp = EditErrorKeyUp
        end
        object EditDump: TEdit
          Left = 152
          Top = 256
          Width = 145
          Height = 21
          TabOrder = 6
          Text = 'dump.dat'
          OnKeyUp = EditDumpKeyUp
        end
        object EditV: TEdit
          Left = 152
          Top = 92
          Width = 145
          Height = 21
          TabOrder = 2
          Text = 'vv.dat'
          OnKeyUp = EditVKeyUp
        end
      end
      object GroupBox2: TGroupBox
        Left = 24
        Top = 24
        Width = 249
        Height = 233
        Caption = 'Ввод начальных данных'
        TabOrder = 1
        object Label8: TLabel
          Left = 32
          Top = 40
          Width = 60
          Height = 13
          Caption = 'Имя файла:'
        end
        object EditInit: TEdit
          Left = 32
          Top = 64
          Width = 97
          Height = 21
          TabOrder = 0
          Text = 'init.dat'
          OnKeyUp = EditInitKeyUp
        end
        object InputInit: TButton
          Left = 32
          Top = 120
          Width = 169
          Height = 41
          Caption = 'Ввести данные в программу'
          TabOrder = 1
          OnClick = InputInitClick
        end
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Calculation'
      object BeginWork: TButton
        Left = 480
        Top = 56
        Width = 145
        Height = 45
        Caption = 'Начать счет'
        ParentShowHint = False
        ShowHint = False
        TabOrder = 1
        OnClick = BeginWorkClick
      end
      object Screen: TMemo
        Left = 24
        Top = 48
        Width = 417
        Height = 329
        Color = clCaptionText
        ScrollBars = ssVertical
        TabOrder = 0
      end
      object DoDump: TButton
        Left = 480
        Top = 227
        Width = 145
        Height = 45
        Caption = 'Сбросить данные в файл'
        TabOrder = 2
        OnClick = DoDumpClick
      end
      object Exit: TButton
        Left = 480
        Top = 312
        Width = 145
        Height = 45
        Caption = 'Окончание работы'
        TabOrder = 3
        OnClick = ExitClick
      end
      object ButtonContinue: TButton
        Left = 480
        Top = 141
        Width = 145
        Height = 45
        Caption = 'Продолжить счет'
        TabOrder = 4
        Visible = False
        OnClick = ContinueWorkClick
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Printing'
      object Label9: TLabel
        Left = 384
        Top = 56
        Width = 184
        Height = 13
        Caption = 'Как часто проверять расходимость:'
      end
      object Label10: TLabel
        Left = 416
        Top = 144
        Width = 40
        Height = 13
        Caption = 'каждый'
      end
      object Label11: TLabel
        Left = 624
        Top = 85
        Width = 81
        Height = 13
        Caption = 'шаг по времени'
      end
      object Label12: TLabel
        Left = 624
        Top = 144
        Width = 81
        Height = 13
        Caption = 'шаг по времени'
      end
      object Label13: TLabel
        Left = 416
        Top = 85
        Width = 40
        Height = 13
        Caption = 'каждый'
      end
      object Label14: TLabel
        Left = 384
        Top = 120
        Width = 153
        Height = 13
        Caption = 'Как часто выводить на экран:'
      end
      object CheckScreenOutput: TCheckBox
        Left = 72
        Top = 88
        Width = 169
        Height = 17
        Hint = 'удобно'
        Caption = 'Выводить данные на экран'
        Checked = True
        ParentShowHint = False
        ShowHint = True
        State = cbChecked
        TabOrder = 1
      end
      object CheckFileOutput: TCheckBox
        Left = 72
        Top = 56
        Width = 177
        Height = 17
        Hint = 'без этого нельзя'
        Caption = 'Выводить данные в файлы'
        Checked = True
        ParentShowHint = False
        ShowHint = True
        State = cbChecked
        TabOrder = 0
      end
      object CheckClear: TCheckBox
        Left = 72
        Top = 120
        Width = 249
        Height = 17
        Caption = 'Очищать экран после вывода одного шага'
        TabOrder = 2
      end
      object CSpinChecking: TCSpinEdit
        Left = 480
        Top = 80
        Width = 121
        Height = 22
        TabStop = True
        MaxValue = 1000
        MinValue = 1
        ParentColor = False
        TabOrder = 3
        Value = 10
        OnChange = CSpinCheckingChange
      end
      object CSpinOutput: TCSpinEdit
        Left = 480
        Top = 139
        Width = 121
        Height = 22
        TabStop = True
        MaxValue = 1000
        MinValue = 1
        ParentColor = False
        TabOrder = 4
        Value = 20
        OnChange = CSpinOutputChange
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Parameters'
      object Label15: TLabel
        Left = 96
        Top = 68
        Width = 20
        Height = 13
        Caption = 'Re='
      end
      object EditRe: TEdit
        Left = 120
        Top = 64
        Width = 121
        Height = 21
        TabOrder = 0
        Text = '10000'
        OnKeyUp = EditReKeyUp
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Tasks'
    end
  end
end
