object MainWindow: TMainWindow
  Left = 1
  Top = 1
  Width = 636
  Height = 478
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
    Width = 798
    Height = 460
    Anchors = [akLeft, akTop, akRight, akBottom]
    PageIndex = 4
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
      object GroupOutput: TGroupBox
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
          OnChange = EditDumpChange
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
      object GroupInit: TGroupBox
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
          Top = 168
          Width = 169
          Height = 41
          Caption = 'Ввести данные в программу'
          TabOrder = 1
          OnClick = InputInitClick
        end
        object RadioInitVel: TRadioButton
          Left = 32
          Top = 104
          Width = 169
          Height = 17
          Caption = 'Ввести начальную скорость'
          Checked = True
          TabOrder = 2
          TabStop = True
        end
        object RadioInitContinue: TRadioButton
          Left = 32
          Top = 136
          Width = 113
          Height = 17
          Caption = 'Продолжить счет'
          TabOrder = 3
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
      object ButtonRecontinue: TButton
        Left = 480
        Top = 141
        Width = 145
        Height = 45
        Caption = 'Начать счет сызнова'
        TabOrder = 4
        Visible = False
        OnClick = RecontinueClick
      end
      object ClearScr: TButton
        Left = 24
        Top = 392
        Width = 65
        Height = 25
        Caption = 'Очистить'
        TabOrder = 5
        OnClick = ClearScrClick
      end
      object SelPrior: TComboBox
        Left = 656
        Top = 152
        Width = 113
        Height = 21
        Hint = 'по умолчанию самый низкий'
        ItemHeight = 13
        ParentShowHint = False
        ShowHint = True
        TabOrder = 6
        Text = 'Приоритет'
        OnChange = SelPriorChange
        Items.Strings = (
          '0.в безделье'
          '1.низкий'
          '2.ниже среднего'
          '3.нормальный'
          '4.выше среднего'
          '5.высокий'
          '6.реального времени')
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
        Width = 208
        Height = 13
        Caption = 'Как часто выводить на экран и в файлы:'
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
        Checked = True
        State = cbChecked
        TabOrder = 2
      end
      object EditCheckStep: TEdit
        Left = 480
        Top = 80
        Width = 120
        Height = 21
        TabOrder = 4
        Text = '10'
        OnChange = EditCheckStepChanging
      end
      object UpDownCheck: TUpDown
        Left = 600
        Top = 80
        Width = 15
        Height = 21
        Associate = EditCheckStep
        Min = 1
        Position = 10
        TabOrder = 3
        Wrap = False
      end
      object EditOutStep: TEdit
        Left = 480
        Top = 144
        Width = 120
        Height = 21
        TabOrder = 5
        Text = '20'
        OnChange = EditOutStepChanging
      end
      object UpDownOut: TUpDown
        Left = 600
        Top = 144
        Width = 15
        Height = 21
        Associate = EditOutStep
        Min = 1
        Max = 1000
        Position = 20
        TabOrder = 6
        Wrap = False
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
      object GroupTurbulent: TGroupBox
        Left = 40
        Top = 264
        Width = 425
        Height = 153
        Caption = 'Турбулентная вязкость'
        TabOrder = 0
        object Label18: TLabel
          Left = 24
          Top = 60
          Width = 180
          Height = 13
          Caption = 'Коэффициент пропорциональности'
        end
        object CheckNut: TCheckBox
          Left = 16
          Top = 32
          Width = 369
          Height = 17
          Caption = 'Использовать турбулентную вязкость через каскадные функции'
          Checked = True
          State = cbChecked
          TabOrder = 0
          OnClick = CheckNutClick
        end
        object EditMaschtab: TEdit
          Left = 224
          Top = 56
          Width = 121
          Height = 21
          TabOrder = 1
          Text = '1000'
          OnChange = EditMaschtabChange
        end
      end
      object EditRe: TEdit
        Left = 120
        Top = 64
        Width = 121
        Height = 21
        TabOrder = 1
        Text = '10000'
        OnKeyUp = EditReKeyUp
      end
      object GroupInitVariant: TGroupBox
        Left = 40
        Top = 104
        Width = 305
        Height = 145
        Caption = 'Выберите тип начальных условий для скорости'
        TabOrder = 2
        object RadioZeroInit: TRadioButton
          Left = 24
          Top = 24
          Width = 113
          Height = 17
          Caption = 'Нулевая'
          TabOrder = 0
        end
        object RadioParabole: TRadioButton
          Left = 24
          Top = 56
          Width = 129
          Height = 17
          Caption = 'Профиль Пуазейля'
          Checked = True
          TabOrder = 1
          TabStop = True
        end
      end
      object GroupNoise: TGroupBox
        Left = 496
        Top = 264
        Width = 217
        Height = 153
        Caption = 'Параметры шума'
        TabOrder = 3
        object Label16: TLabel
          Left = 24
          Top = 56
          Width = 78
          Height = 13
          Caption = 'Величина шума'
        end
        object Label17: TLabel
          Left = 24
          Top = 116
          Width = 80
          Height = 13
          Caption = 'Величина возм.'
        end
        object CheckNoise: TCheckBox
          Left = 24
          Top = 24
          Width = 97
          Height = 17
          Caption = 'Включить шум'
          Checked = True
          State = cbChecked
          TabOrder = 0
          OnClick = CheckNoiseClick
        end
        object EditNoise: TEdit
          Left = 120
          Top = 52
          Width = 65
          Height = 21
          TabOrder = 1
          Text = '0.1'
          OnKeyUp = EditNoiseKeyUp
        end
        object CheckNorm: TCheckBox
          Left = 24
          Top = 88
          Width = 177
          Height = 17
          Caption = 'Включить норм. возмущения'
          TabOrder = 2
          OnClick = CheckNormClick
        end
        object EditNorm: TEdit
          Left = 120
          Top = 112
          Width = 65
          Height = 21
          Enabled = False
          TabOrder = 3
          Text = '0'
          OnChange = EditNormChange
        end
      end
      object GroupGrid: TGroupBox
        Left = 392
        Top = 104
        Width = 353
        Height = 145
        Caption = 'Параметры сетки'
        TabOrder = 4
        object Label19: TLabel
          Left = 16
          Top = 24
          Width = 89
          Height = 13
          Caption = 'Кол-во узлов по x'
        end
        object Label20: TLabel
          Left = 16
          Top = 100
          Width = 89
          Height = 13
          Caption = 'Кол-во узлов по z'
        end
        object Label21: TLabel
          Left = 16
          Top = 62
          Width = 89
          Height = 13
          Caption = 'Кол-во узлов по y'
        end
        object EditNx: TEdit
          Left = 144
          Top = 20
          Width = 65
          Height = 21
          TabOrder = 0
          Text = 'EditNx'
        end
        object EditNy: TEdit
          Left = 144
          Top = 59
          Width = 65
          Height = 21
          TabOrder = 1
          Text = 'EditNy'
        end
        object EditNz: TEdit
          Left = 144
          Top = 96
          Width = 65
          Height = 21
          TabOrder = 2
          Text = 'EditNz'
        end
        object ButtonSaveNi: TButton
          Left = 224
          Top = 45
          Width = 113
          Height = 49
          Caption = 'Принять изменения'
          TabOrder = 3
          OnClick = ButtonSaveNiClick
        end
      end
    end
    object TTabPage
      Left = 4
      Top = 24
      Caption = 'Tasks'
    end
  end
  object StatusBarCalc: TStatusBar
    Left = 0
    Top = 460
    Width = 790
    Height = 20
    Panels = <
      item
        Text = 'Начало работы (нажмите кнопку)'
        Width = 200
      end
      item
        Text = 'Ничего не сделано'
        Width = 50
      end>
    SimplePanel = False
  end
end
